from collections import defaultdict, namedtuple
import os
from uuid import uuid4
import re
import time

from flask import (
    Blueprint, flash, request, redirect, render_template, url_for, current_app, send_file, make_response, Response)
from flask import current_app as app
from glob import glob
from werkzeug.utils import secure_filename # to secure file
from .worker import run_task
import io
import json
import csv
import pathogenprofiler as pp
from .files import get_paired_fastq_samples, get_single_fasta_samples, get_single_fastq_samples
import shutil

bp = Blueprint('main', __name__)

def get_upload_dir(upload_id):
    return os.path.join(app.config["UPLOAD_FOLDER"],upload_id)

@bp.route('/')
def index():
    return render_template("pages/index.html")

@bp.route('/upload', methods=["GET", "POST"])
def upload():
    random_id = str(uuid4())
    # random_id = '9b5adba5-4602-428e-ac48-4bef57948028'
    if request.method == "POST":
        platform = request.form["radio_platform"]
        filetype = request.form['radio_filetype']
        forward_regex = "(.+)"+request.form['forward-suffix'] if request.form['forward-suffix']!="" else '(.+)_1.fastq.gz'
        reverse_regex = "(.+)"+request.form['reverse-suffix'] if request.form['reverse-suffix']!="" else '(.+)_2.fastq.gz'

        print(filetype)
        print(request.form)
        runs = []
        upload_id = request.form['submit_button']
        upload_dir = get_upload_dir(upload_id)
        if not os.path.isdir(upload_dir):
            flash("No new files uploaded","danger")
            return render_template("pages/upload.html", random_id=random_id)

        new_upload_id = str(uuid4())
        new_upload_dir = get_upload_dir(new_upload_id)
        os.rename(upload_dir,new_upload_dir)

        files = glob(f'{new_upload_dir}/*')
        print(files)
        if filetype=='paired':
            samples = get_paired_fastq_samples(files,r1_suffix=forward_regex,r2_suffix=reverse_regex)
        elif filetype=='single':
            samples = get_single_fastq_samples(files,r1_suffix='(.+).fastq.gz')
        elif filetype=='fasta':
            samples = get_single_fasta_samples(files,fasta_suffix='(.+).fasta')
        print(samples)
        if len(samples)==0:
            flash("No valid files found. Check and see if your file suffix is correct.","danger")
            return render_template("pages/upload.html", random_id=random_id)
        for s in samples:
            print("*"*100)
            print(s)

            run_id = str(uuid4())
            log_file = "%s/%s.log" % (app.config["RESULTS_DIR"], run_id) 
            print("Logging to:",log_file)
            with open(log_file, "w") as O:
                O.write("Submitting job: %s\n" % run_id)

            run_task.delay(
                files=s.files, 
                filetype=s.filetype, 
                platform=platform,
                run_id=run_id, 
                results_dir=app.config["RESULTS_DIR"],
                threads=app.config["THREADS"]
            )

            runs.append({"id":run_id, "files":s.files})
        analysis_id = str(uuid4())
        with open("%s/%s.json" % (app.config["RESULTS_DIR"],analysis_id), "w") as O:
            json.dump(runs,O)
        return redirect(url_for("main.upload_runs_id", analysis_id=analysis_id))
        
    
    return render_template("pages/upload.html", random_id=random_id)

file_patterns = {
    "fasta": "\.fasta$|\.fa$|\.fna$",
    "fastq": "\.fastq.[A-Za-z]*$|\.fq.[A-Za-z]*$",
}

@bp.route('/run_result/<uuid:analysis_id>')
def upload_runs_id(analysis_id):
    data = json.load(open("%s/%s.json" % (app.config["RESULTS_DIR"], analysis_id)))
    print(data)
    for d in data:
        d["link"] = '<a href="' + url_for("main.result_id", run_id=d["id"]) + '">' + d["id"] + '</a>'
        d["files"] = ", ".join([x.split("/")[-1] for x in d["files"]])
    return render_template('/pages/analysis_id.html', runs = data)

def get_filetype(filename):
    for key,pattern in file_patterns.items():
        if re.search(pattern,filename.strip().lower()):
            print(key)
            return key
    return None

# def get_files_in_dir(upload_dir):
#     file_list = glob("%s/*" % upload_dir)
#     File = namedtuple("File", "files type")
#     files = []
#     fastqs = []
#     for f in file_list:
#         ftype = get_filetype(f)
#         if ftype=="fastq":
#             fastqs.append(f)
#         else:
#             files.append(File((f,),ftype))

#     pattern = "(.*)(_R?[12])(\.fastq.[A-Za-z]*$)|(.*)(_R?[12])(\.fq.[A-Za-z]*$)"

#     fastqs = sorted(fastqs)
#     while len(fastqs)>0:
#         f = fastqs.pop(0)
#         r = re.search(pattern,f)
#         if r:
#             if r.group(2)=="_1":
#                 potential_pair = r.group(1)+"_2"+r.group(3)
#             else:
#                 potential_pair = r.group(1)+"_R2"+r.group(3)

#             print(f"Looking for {potential_pair} in {str(fastqs)}: {potential_pair in fastqs}")
#             if potential_pair in fastqs:
#                 pair = fastqs.pop(fastqs.index(potential_pair))
#                 files.append(File((f,pair),"fastq"))
#             else:
#                 files.append(File((f,),"fastq"))
#         else:
#             files.append(File((f,),"fastq"))

#     return(files)

def is_legal_filetype(filename):
    if get_filetype(filename):
        return True
    else:
        return False

def get_conf(results):
    db_name = results['resistance_db']['name']
    conf = pp.get_db('malaria_profiler',db_name)
    return conf

def parse_result_summary(json_file):
    data = {}
    results = json.load(open(json_file))
    data['species'] = ", ".join([e['species'] for e in results['species']['species']])
    if 'barcode' in results:
        data['subspecies'] = ", ".join([e['id'] for e in results['barcode']])
    return data

def add_linebreaks(text):
    # make each line a new div
    return "<div>"+"</div><div>".join(text.split("\n")) + "</div>"

def get_drug_table(dr_variants,dr_genes,conf=None,drugs=None):
    all_drugs = conf['drugs']
    new_table = []
    for v in dr_variants:
        for d in v['drugs']:
            new_row = {
                'drug': d['drug'],
                'gene': v['gene_name'],
                'change': v['change'],
            }
            new_table.append(new_row)
    for g in dr_genes:
        for d in g['drugs']:
            new_row = {
                'drug': d['drug'],
                'gene': g['gene_name'],
                'change': g['type'],
            }
            new_table.append(new_row)
    variant_drugs = list(set([r['drug'] for r in new_table]))
    for d in all_drugs:
        if d not in variant_drugs:
            new_table.append({
                'drug': d,
                'gene': '',
                'change': '',
                'confidence': '',
                'comment': '',
            })
    
    new_table = [r for r in new_table if r['drug'] in all_drugs]
    new_table = sorted(new_table, key=lambda x: all_drugs.index(x['drug']))
    for drug in all_drugs:
        drugrows = [d for d in new_table if d['drug'] == drug]
        for i,r in enumerate(drugrows):
            if i == 0:
                r['drug-rowspan'] = len(drugrows)
            generows = [g for g in drugrows if g['gene'] == r['gene']]
            for j,g in enumerate(generows):
                if j == 0:
                    g['gene-rowspan'] = len(generows)
    return new_table

def get_reference_files(conf):
    local_ref_file_name = url_for('static', filename='reference_files/%s' % conf['ref'].split("/")[-1])
    if not os.path.isfile(local_ref_file_name):
        # copy the file to the static folder
        shutil.copy(conf['ref'], os.path.join(app.config["REFERENCE_DIR"],conf['ref'].split("/")[-1]))
    
    local_ref_file_index_name = url_for('static', filename='reference_files/%s.fai' % conf['ref'].split("/")[-1])
    if not os.path.isfile(local_ref_file_index_name):
        # copy the file to the static folder
        shutil.copy(conf['ref']+".fai", os.path.join(app.config["REFERENCE_DIR"],conf['ref'].split("/")[-1]+".fai"))
    
    local_gff_file_name = url_for('static', filename='reference_files/%s' % conf['gff'].split("/")[-1])
    if not os.path.isfile(local_gff_file_name):
        # copy the file to the static folder
        shutil.copy(conf['gff'], os.path.join(app.config["REFERENCE_DIR"],conf['gff'].split("/")[-1]))
    return {
        'ref': local_ref_file_name,
        'gff': local_gff_file_name
    }

@bp.route('/result/<uuid:run_id>')
def result_id(run_id):
    log_file = "%s/%s.log" % (app.config["RESULTS_DIR"], run_id)
    if not os.path.isfile(log_file):
        flash("Error! Result with ID:%s doesn't exist" % run_id, "danger")
        return render_template('pages/result.html')
    json_file = "%s/%s.results.json" % (app.config["RESULTS_DIR"], run_id)
    log_text = add_linebreaks(open(log_file).read())

    if not os.path.isfile(json_file):
        
        return render_template('pages/still-waiting.html', run_id=run_id, log_text = log_text)
    else:
        results = json.load(open(json_file))
        data = parse_result_summary(json_file)
        if results['result_type']=='Species':
            return render_template('pages/species-result.html', run_id=run_id, results = results, data = data, log_text=log_text)
        else:
            db_name = results['resistance_db']['name']
            conf = pp.get_db('ntm-profiler',db_name)
            reference = get_reference_files(conf)
            data['drug_resistance_table'] = get_drug_table(results['dr_variants'],results['dr_genes'],conf)
            files = {}
            bam_filename = "%s/%s.bam" % (app.config["RESULTS_DIR"], run_id)
            if os.path.isfile(bam_filename):
                files['bam'] = bam_filename

            return render_template('pages/full-result.html', run_id=run_id, results=results, data = data, log_text=log_text,reference=reference,files=files)

@bp.route('/result/<uuid:run_id>/download', methods=['GET', 'POST'])
def download(run_id):
        result_file = "%s/%s.results.txt" % (app.config["RESULTS_DIR"], run_id)
        return send_file(result_file, as_attachment=True)

@bp.route('/result', methods=['POST', 'GET'])
def result():
    if request.method == "POST":
        if "result_submit" in request.form:
            run_id = request.form["result_id"].strip()
            return redirect(url_for('main.result_id', run_id=run_id))
    return render_template("pages/result.html")



@bp.route('/file_upload/<uuid:upload_id>',methods=('GET','POST'))
def file_upload(upload_id):

    upload_id = str(upload_id)
    file = request.files['file']
    if not is_legal_filetype(file.filename):
        return make_response(('Unknown file type', 400))

    upload_dir = os.path.join(app.config["UPLOAD_FOLDER"],upload_id)
    if not os.path.isdir(upload_dir):
        os.mkdir(upload_dir)
    save_path = os.path.join(upload_dir, file.filename)
    current_chunk = int(request.form['dzchunkindex'])
    # If the file already exists it's ok if we are appending to it,
    # but not if it's new file that would overwrite the existing one
    if os.path.exists(save_path) and current_chunk == 0:
        # 400 and 500s will tell dropzone that an error occurred and show an error
        return make_response(('File already exists', 400))
    try:
        with open(save_path, 'ab') as f:
            f.seek(int(request.form['dzchunkbyteoffset']))
            f.write(file.stream.read())
    except OSError:
        # log.exception will include the traceback so we can see what's wrong
        # log.exception('Could not write to file')
        return make_response(("Not sure why,"
                              " but we couldn't write the file to disk", 500))
    total_chunks = int(request.form['dztotalchunkcount'])
    if current_chunk + 1 == total_chunks:
        # This was the last chunk, the file should be complete and the size we expect
        if os.path.getsize(save_path) != int(request.form['dztotalfilesize']):
            return make_response(('Size mismatch', 500))
        else:
            print(f'File {file.filename} has been uploaded successfully from session {upload_id} to {save_path}')
    else:
        print(f'Chunk {current_chunk + 1} of {total_chunks} for file {file.filename} complete')
    return make_response(("Chunk upload successful", 200))