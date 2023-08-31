# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz


configfile: "config/config.yaml"
BASE = "/agr/persist/projects/2023-bjp-rmre-seq/SMK-rMRE-Seq"

LIBRARY = config["library"]

# def get_library_barcodes(wildcards) # TODO Input function that accepts library, and returns the barcodes to expand with
BARCODES = ['BC001','BC002','BC003','BC004','BC005','BC006','BC007','BC008','BC009','BC010','BC011','BC012','BC013','BC014','BC015','BC016','BC017','BC018','BC019','BC020','BC021','BC022','BC023','BC024']

import os
import pandas as pd
from glob import glob


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


rule all:
    input:
        expand("results/{library}/02_align/{library}.{barcode}.sam", library = LIBRARY, barcode = BARCODES),
        # expand("results/{library}/00_stats/{library}.{barcode}.samtools.stats.txt", library = LIBRARY, barcode = BARCODES),


rule bowtie2:
    input:
        demuxed = "results/{library}/01_cutadapt/{library}.{barcode}.fastq.gz"
    output:
        aligned = "results/{library}/02_align/{library}.{barcode}.sam",
    log:
        "logs/bowtie2/bowtie2.{library}.{barcode}.log"
    benchmark:
        'benchmarks/bowtie2.{library}.{barcode}.txt'
    conda:
        'bowtie2-2.5.1'
    threads: 16
    params:
        bowtie2_reference = "resources/reference/GCF_016772045.1/GCF_016772045.1_ARS-UI_Ramb_v2.0", #TODO pass via config for CLI integration
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 12),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition="compute"
    shell:
        "bowtie2 "
        "--seed 1953 "
        "--threads {threads} "
        "-x {params.bowtie2_reference} "
        "-U {input.demuxed} "
        "-S {output.aligned} "
        " | tee {log} "


rule bam_stats:
