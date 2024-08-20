import os
from celery import shared_task
from celery.result import AsyncResult
import json

import subprocess as sp #TODO check multiprocessing
from celery.utils.log import get_task_logger #TODO logger
import time 
import shutil
import pathogenprofiler as pp
import sys
from typing import List
import logging

BASH_TIMEOUT = os.environ.get('BASH_TIMEOUT', 1200)






def get_status(run_id):
    return AsyncResult(run_id).state

@shared_task
def run_task(
        files: List[str], 
        filetype: str,
        platform: str, 
        run_id: str, 
        results_dir: str, 
        threads: int = 1
    ):
    # log to file
    log_fh = open(f"{results_dir}/{run_id}.log")

    print(f"Writing log to: {results_dir}/{run_id}.log")


    print(filetype)
    if filetype=='fasta':
        cmd = f"ntm-profiler profile -f {files[0]} --dir {results_dir} --prefix {run_id} --platform {platform} -t {threads} --txt " 
    elif filetype=='paired-fastq':
        cmd = f"ntm-profiler profile -1 {files[0]} -2 {files[1]} --dir {results_dir} --prefix {run_id} --platform {platform} -t {threads} --txt " 
    elif filetype=='single-fastq':
        cmd = f"ntm-profiler profile -1 {files[0]} --dir {results_dir} --prefix {run_id} --platform {platform} -t {threads} --txt " 
    cmd += f" | tee {results_dir}/{run_id}.log"
    print(cmd)

    sp.run(cmd, shell=True)
    if os.path.isfile(f"{results_dir}/{run_id}.bam"):
        db_name = json.load(open(f"{results_dir}/{run_id}.results.json"))["resistance_db"]["name"]
        bed_file = f"{sys.base_prefix}/share/ntm-profiler/{db_name}.bed"
        print(bed_file)
        sp.call(f"samtools view -bL {bed_file} {results_dir}/{run_id}.bam > {results_dir}/{run_id}.bed.bam", shell=True)
        sp.call(f"mv {results_dir}/{run_id}.bed.bam {results_dir}/{run_id}.bam", shell=True)
        sp.call(f"samtools index {results_dir}/{run_id}.bam", shell=True)

    

@shared_task
def remote_profile(ftype, files, run_id, results_dir, platform, species, threads = 1):
    tmp_dir = f"/tmp/runs/"
    conf = {
        "run_id": run_id,
        "ftype": ftype,
        "platform": platform,
        "files": files,
        "species": species
    }
    run_file = f"{tmp_dir}/{run_id}.run_file.json"
    json.dump(conf,open(run_file,"w"))
    server_result_file = f"{tmp_dir}/{run_id}.completed.json"
    while True:
        time.sleep(1)
        if os.path.exists(server_result_file):
            break
    server_result = json.load(open(server_result_file,"r"))
    print(server_result)
    for val in server_result.values():
        local_file_name = val.split("/")[-1]
        shutil.copyfile(f"{tmp_dir}/{local_file_name}",f"{results_dir}/{local_file_name}" )
        os.remove(f"{tmp_dir}/{local_file_name}")
    sp.call(f"samtools index {results_dir}/{run_id}.bam", shell=True)
    os.remove(server_result_file)
    os.remove(run_file)