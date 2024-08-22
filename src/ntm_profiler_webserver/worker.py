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
    
    # remove the files
    for f in files:
        os.remove(f)

    with open(f"{results_dir}/{run_id}.log","a") as LOG:
        LOG.write("\nDONE\n")
