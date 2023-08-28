# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"

import os
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


rule cutadapt: # demultiplexing GBS reads
    input:
        barcodes = "resources/barcodes.fasta",
        library = "data/{library}.fastq.gz",
    output:
        demuxed = directory("results/{library}/01_cutadapt"),
    log:
        "logs/cutadapt.{library}.log.json"
    container:
        "docker://quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1"
    benchmark:
        "benchmarks/cutadapt.{library}.txt"
    threads: 32
    resources:
        mem_gb=8,
        time="01:00:00",
	partition="compute"
    shell:
        "mkdir -p {output.demuxed} && "
        "zcat {input.library} | "
        "cutadapt "
        "-json={log} "
        "-j {threads} "
        "--discard-untrimmed "
        "--length 65 "
        "-m 50 "
        "--action retain "
        "-e 0 "
        "--no-indels "
        "-g ^file:{input.barcodes} "
        r"-o '{output.demuxed}/{{name}}.fastq.gz' "
        "- && "
        " exit 0;"

