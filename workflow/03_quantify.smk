# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: 'config/config.yaml'

import os
import pandas as pd


onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')
    print("Found: ")


wildcard_constraints:
    samples="\w+"


FIDs, = glob_wildcards('results/01_cutadapt/{samples}.fastq.gz')


rule all:
    input:
        'results/00_QC/seqkit.report.raw.txt',
        'results/00_QC/seqkit.report.prepared.txt',


checkpoint seqkitRaw:
    input:
        expand('results/01_prepare/{samples}.sana.fastq.gz', samples = FIDs),
    output:
        'results/00_QC/seqkit.report.raw.txt'
    benchmark:
        'benchmarks/seqkitRaw.txt'
    #container:
    #    'docker://quay.io/biocontainers/seqkit:2.2.0--h9ee0642_0' 
    conda:
        #'env/seqkit.yaml'
        'seqkit'
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 4 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) * 60),
        partition="large,milan"
    shell:
        'seqkit stats -j {threads} -a {input} > {output} '


# STANDARD READ FILTERING AND QC RULES
rule bbduk:
    input:
        reads = 'results/01_prepare/{samples}.sana.fastq.gz',
    output:
        bbdukReads = temp('results/01_prepare/{samples}.bbduk.fastq.gz')
    log:
        'logs/bbduk/bbduk.{samples}.log'
    benchmark:
        'benchmarks/bbduk.{samples}.txt'
    conda:
        'bbduk'
    threads: 8
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) * 4),
        time = lambda wildcards, attempt: 8 + ((attempt - 1) * 10),
        partition='compute',
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.reads} '
        'entropy=0.3 '
        'entropywindow=50 '
        'trimpolygright=5 '
        'qtrim=r '
        'trimq=20 '
        'out={output.bbdukReads} '
        '2>&1 | tee {log}'

