# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"

import os
from glob import glob

BASE = "/agr/persist/projects/2023-bjp-rmre-seq/SMK-rMRE-Seq"

LIBRARY = config["library"]



#SAMPLES = function_parse_keyfile() #TODO


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


rule all:
    input:
        expand("results/{library}/01_cutadapt", library = LIBRARY),
        expand("results/{library}/00_fastqc", library = LIBRARY),



rule optical_duplicates: # Removing optical duplicates from NovaSeq run
    input:
        reads = "data/{library}.fastq.gz",
    output:
        dedupe = "results/{library}/00_dedupe/{library}.clumpify-15000.fastq.gz", #TODO Find a way to link with params; SMK issue on Github requesting this currently
    log:
        "logs/clumpify/{library}.clumpify.log"
    benchmark:
        "benchmarks/clumpify/clumpify.{library}.txt"
    conda:
        "bbmap-39.01"
    threads: 24
    resources:
        mem_gb = lambda wildcards, attempt: 80 + ((attempt - 1) * 80),
        time = lambda wildcards, attempt: 240 + ((attempt - 1) * 120),
	    partition="compute"
    params:
        distance = "15000" #TODO Move to config
    shell:
        "clumpify.sh "
        "-Xmx{resources.mem_gb}g "
        "optical "
        "dedupe "
        "dupedist={params.distance} "
        "subs=0 "
        "in={input.reads} "
        "out={output.dedupe} "
        "2> {log} "


rule adapter_trimming: # Removing Illumina sequencing adapters
    input:
        dedupe = "results/{library}/00_dedupe/{library}.clumpify-15000.fastq.gz", 
    output:
        trimmed = "results/{library}/00_trimmed/{library}.clumpify-15000.trimmed.fastq.gz",
    log:
        "logs/bbduk/{library}.bbduk.trimming.log"
    conda:
        "bbmap-39.01"
    threads: 24
    resources:
        mem_gb = lambda wildcards, attempt: 80 + ((attempt - 1) * 80),
        time = lambda wildcards, attempt: 240 + ((attempt - 1) * 120),
	    partition="compute"
    params:
        adapters = "resources/illumina_adapters.fasta"
    shell:
        "bbduk.sh "
        "-Xmx{resources.mem_gb}g "
        "stats={log} "
        "in={input.dedupe} "
        "adapters={params.adapters} "
        "ktrim=r "
        "ktrim=r "
        "k=23 "
        "mink=8 "
        "hdist=1 "
        "trimpolygright=9 "
        "out={output.trimmed} "


rule cutadapt: # demultiplexing GBS reads
    input:
        barcodes = "resources/barcodes.fasta",
        trimmed = "results/{library}/00_trimmed/{library}.clumpify-15000.trimmed.fastq.gz",
    output:
        demuxed = directory("results/{library}/01_cutadapt"),
    log:
        "logs/cutadapt/{library}.cutadapt.log"
    benchmark:
        "benchmarks/cutadapt.{library}.txt"
    conda:
        "cutadapt-4.4"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 240 + ((attempt - 1) * 120),
	    partition="compute"
    shell:
        "mkdir -p {output.demuxed} && "
        "zcat {input.trimmed} | "
        "cutadapt "
        #"--json={log} "
        "-j {threads} "
        "--discard-untrimmed "
        "--length 65 "
        "-m 25 "
        #"--action=retain "
        "-e 0 "
        "-O 5 "
        "--no-indels "
        "-g ^file:{input.barcodes} "
        r'-o "{output.demuxed}/{wildcards.library}.{{name}}.fastq.gz" '
        "- "
        " | tee {log} "


rule fastqc_reads:
    input:
        demuxed = "results/{library}/01_cutadapt",
    output:
        fastqc = directory("results/{library}/00_fastqc")
    log:
        "logs/fastqc/{library}.fastqc.log"
    conda:
        "fastqc-0.12.1"
    threads:
        24
    resources:
        mem_gb = lambda wildcards, attempt: 24 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 240 + ((attempt - 1) * 120),
	    partition="compute"
    shell:
        "mkdir -p {output.fastqc} && "
        "fastqc "
        "-t 24 "
        "-o {output.fastqc} "
        "{input.demuxed}/*.fastq.gz "
        " | tee {log} "


