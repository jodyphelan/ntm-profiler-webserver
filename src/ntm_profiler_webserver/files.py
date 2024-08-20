import os
import subprocess as sp
import re
from collections import defaultdict
from typing import List

__version__ = "0.0.2"

class Sample:
    def __init__(self,prefix: str,r1: str) -> None:
        self.prefix = prefix


class FastqSingleSample(Sample):
    """
    Class to represent a sample with a single fastq file
    
    Parameters
    ----------
    prefix : str
        Sample ID
    r1 : str
        Path to the R1 file
    
    Attributes
    ----------
    prefix : str
        Sample ID
    r1 : str
        Path to the R1 file
    """
    filetype='single-fastq'
    def __init__(self,prefix: str,r1: str) -> None:
        self.prefix = prefix
        self.files = [r1]

    def __repr__(self) -> str:
        return f"Sample(prefix={self.prefix},r1={self.files[0]})"

class FastaSample(Sample):
    """
    Class to represent a sample with a single fastq file
    
    Parameters
    ----------
    prefix : str
        Sample ID
    r1 : str
        Path to the R1 file
    
    Attributes
    ----------
    prefix : str
        Sample ID
    r1 : str
        Path to the R1 file
    """
    filetype='fasta'
    def __init__(self,prefix: str,fasta: str) -> None:
        self.prefix = prefix
        self.files = [fasta]

    def __repr__(self) -> str:
        return f"Sample(prefix={self.prefix},fasta={self.files[0]})"


class FastqPairedSample(Sample):
    """
    Class to represent a sample with R1 and R2 files
    
    Parameters
    ----------
    prefix : str
        Sample ID
    r1 : List[str]
        List of R1 files
    r2 : List[str]
        List of R2 files
    
    Attributes
    ----------
    prefix : str
        Sample ID
    r1 : List[str]
        List of R1 files
    r2 : List[str]
        List of R2 files
    multiple : bool
        Whether there are multiple R1 and R2 files
    """
    filetype='paired-fastq'
    def __init__(self,prefix: str,r1: str,r2: str) -> None:
        self.prefix = prefix
        self.files = [r1,r2]

    def __repr__(self) -> str:
        return f"Sample(prefix={self.prefix},r1={self.files[0]},r2={self.files[1]})"


def get_single_fastq_samples(files: List[str], r1_suffix: str) -> List[FastqSingleSample]:
    """
    Function to sort out single files from a list of files

    Parameters
    ----------
    files : List[str]
        List of files to sort out
    r1_suffix : str
        Suffix for R1 files

    Returns
    -------
    List[Sample]
        List of Sample objects
    """
    prefixes = defaultdict(list)

    for f in files:
        tmp = re.search("%s" % r1_suffix,f)
        p = None
        if tmp:
            p = tmp.group(1).split("/")[-1]
            prefixes[p].append(f)

    runs = []
    for p,vals in prefixes.items():
        vals.sort()
        
        runs.append(
            FastqSingleSample(p,vals[0])
        )

    return runs


def get_paired_fastq_samples(files: List[str], r1_suffix: str,r2_suffix: str) -> List[FastqPairedSample]:
    """
    Function to sort out paired files from a list of files

    Parameters
    ----------
    files : List[str]
        List of files to sort out
    r1_suffix : str
        Suffix for R1 files
    r2_suffix : str
        Suffix for R2 files

    Returns
    -------
    List[Sample]
        List of Sample objects
    """
    prefixes = defaultdict(lambda:{"r1":[],"r2":[]})

    for f in files:
        tmp1 = re.search("%s" % r1_suffix,f)
        tmp2 = re.search("%s" % r2_suffix,f)
        p = None
        if tmp1:
            p = tmp1.group(1).split("/")[-1]
            prefixes[p]['r1'].append(f)
        elif tmp2:
            p = tmp2.group(1).split("/")[-1]
            prefixes[p]['r2'].append(f)

    runs = []
    for p,vals in prefixes.items():

        if len(vals['r1'])!=len(vals['r2']):
            raise ValueError(f"Number of R1 and R2 files for sample {p} do not match")
        vals['r1'].sort()
        vals['r2'].sort()
        runs.append(
            FastqPairedSample(p,vals['r1'][0],vals['r2'][0])
        )
    return runs


def get_single_fasta_samples(files: List[str], fasta_suffix: str) -> List[FastqSingleSample]:
    """
    Function to sort out single files from a list of files

    Parameters
    ----------
    files : List[str]
        List of files to sort out
    r1_suffix : str
        Suffix for R1 files

    Returns
    -------
    List[Sample]
        List of Sample objects
    """
    prefixes = defaultdict(list)

    print("Looking for files with suffix %s in:" % fasta_suffix)
    print(files)

    for f in files:
        tmp = re.search("%s" % fasta_suffix,f)
        p = None
        if tmp:
            p = tmp.group(1).split("/")[-1]
            prefixes[p].append(f)

    runs = []
    for p,vals in prefixes.items():
        vals.sort()
        
        runs.append(
            FastaSample(p,vals[0])
        )

    return runs


def find_samples(files: List[str], single_suffix: str, paired_r1_suffix: str, paired_r2_suffix: str, fasta_suffix: str) -> List[FastqPairedSample]:
    """
    Function to find samples from a list of files
    
    Parameters
    ----------
    files : List[str]
        List of files to sort out
    single_suffix : str
        Suffix for single files
    paired_r1_suffix : str
        Suffix for R1 files
    paired_r2_suffix : str
        Suffix for R2 files
    fasta_suffix : str
        Suffix for fasta files

    Returns
    -------
    List[Sample]
        List of Sample objects
    """
    samples = []
    samples.extend(get_single_fastq_samples(files,single_suffix))
    samples.extend(get_paired_fastq_samples(files,paired_r1_suffix,paired_r2_suffix))
    samples.extend(get_single_fasta_samples(files,fasta_suffix))

    return samples