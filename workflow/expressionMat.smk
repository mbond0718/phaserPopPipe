#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os

##### Load config and sample sheets #####
configfile: "/gpfs/commons/groups/lappalainen_lab/mbond/work/phaserPopPipe/config/config.yaml"

## Read in samplesheet
samples = pd.read_table(config["samplesheet_expr"])

## Convert all columns to strings
samples = samples.astype(str)

## Unique tissues to run pipeline on
tissues = samples['Tissue'].unique().tolist()

## Define actions on success
onsuccess:
    ## Success message
    print("phaserPopPipe_Expression workflow completed successfully!")

##### Define rules #####

rule all:
    input:
        expand("input/{TISSUE}/ref/{TISSUE}_genes.bed", TISSUE = tissues),
        expand("output/expression/{TISSUE}/{TISSUE}expressionMatrix.txt.bed.gz", TISSUE = tissues),
        expand("output/expression/{TISSUE}/{TISSUE}expressionMatrix.txt.bed.gz.tbi", TISSUE = tissues),
        expand("output/expression/{TISSUE}/{TISSUE}expressionMatrix.txt.gw_phased.bed.gz", TISSUE = tissues),
        expand("output/expression/{TISSUE}/{TISSUE}expressionMatrix.txt.gw_phased.bed.gz.tbi", TISSUE = tissues)

        
rule prepareExpressionMat:
   output:
       ref = "input/{TISSUE}/ref/{TISSUE}_genes.bed"
   params: 
       tissue = lambda wildcards: wildcards.TISSUE,
       bed = config['ref'],
       samplesheet = config['samplesheet'],
       ase_dir = config['aseDir'],
       outpath = config['outpath']
   log:
       err = "input/logs/{TISSUE}_prepMatrix.err",
       out = "input/logs/{TISSUE}_prepMatrix.out"
   resources: 
       mem = 16000
   shell:
       """
        mkdir -p input/{wildcards.TISSUE}/ref
        mkdir -p input/{wildcards.TISSUE}/ae
        module load R
        Rscript --vanilla scripts/prepareExpressionMatrix.R {params.tissue} {params.samplesheet} {params.ase_dir} {params.bed} {params.outpath}
        """
rule phaserExpMat:
   input:
        features = "input/{TISSUE}/ref/{TISSUE}_genes.bed"
   output:
       expressionMat = "output/expression/{TISSUE}/{TISSUE}_expressionMatrix.txt.bed.gz",
       expressionMat = "output/expression/{TISSUE}/{TISSUE}_expressionMatrix.txt.bed.gz.tbi",
       expressionMat = "output/expression/{TISSUE}/{TISSUE}_expressionMatrix.txt.gw_phased.bed.gz",
       expressionMat = "output/expression/{TISSUE}/{TISSUE}_expressionMatrix.txt.gw_phased.bed.gz.tbi"
   params: 
       ae_dir = lambda wildcards: f"input/{wildcards.TISSUE}/ae"
   log:
       err = "output/logs/{TISSUE}_phaserExpMat.err",
       out = "output/logs/{TISSUE}_phaserExpMat.out"
   shell:
       """
        module load BCFtools
        python3 scripts/phaser_expr_matrix_ed.py --gene_ae_dir {params.ae_dir} --features {input.features} --o {output.expressionMat}
        """
