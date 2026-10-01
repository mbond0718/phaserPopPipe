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
pairs = pd.read_table(config["samplesheet_cisvar"])

## Convert all columns to strings
samples = samples.astype(str)

## Unique tissues to run pipeline on
tissues = samples['Tissue'].unique().tolist()

## Tissues and the samples that belong to them: 
tissue_samples = samples.groupby("Tissue")["Sample"].apply(list).to_dict()

## Pairs
pair_list = pairs['Pairs'].unique().tolist()

## Tissues and the variant pairs that are being intersected with 
tissue_pair_list = pairs.to_dict(orient="records")

## Define actions on success
onsuccess:
    ## Success message
    print("phaserPopPipe_CisVar workflow completed successfully!")

##### Define rules #####

rule all:
    input:
        expand("input/{TISSUE}/vcf/{TISSUE}_vcf.gz", TISSUE=tissues),
        expand("input/{TISSUE}/vcf/{TISSUE}_vcf.gz.tbi", TISSUE=tissues),
        expand("input/{TISSUE}/vcf/{TISSUE}_samples.txt", TISSUE=tissues),
        expand("input/{TISSUE}/map/{TISSUE}_samplemap.txt", TISSUE=tissues),
        expand("output/cisvar/{TISSUE}/{TISSUE}_{PAIRS}_pairs.txt",
               zip,
               TISSUE=[x['Tissue'] for x in tissue_pair_list],
               PAIRS=[x['Pairs'] for x in tissue_pair_list])

rule makeVCF:
   input:
        vcf = lambda wildcards: f"/gpfs/commons/groups/lappalainen_lab/mbond/work/phaserPop/input/vcf/{wildcards.TISSUE}/{wildcards.TISSUE}_merged.vcf.gz"  
   output:
       vcf = "input/{TISSUE}/vcf/{TISSUE}_vcf.gz",
       indexedvcf = "input/{TISSUE}/vcf/{TISSUE}_vcf.gz.tbi"
   params: 
       samplesheet = config['samplesheet'],
       samples = lambda wildcards: ",".join(tissue_samples[wildcards.TISSUE])
   log:
       err = "input/logs/{TISSUE}_makeVCF.err",
       out = "input/logs/{TISSUE}_makeVCF.out"
   resources: 
       mem = 16000
   shell:
       """
        mkdir -p input/{wildcards.TISSUE}/vcf
        module load BCFtools
        bcftools view -s {params.samples} -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """
rule sampleMap:
   input:
        vcf = "input/{TISSUE}/vcf/{TISSUE}_vcf.gz"
   output:
       vcfSamples = "input/{TISSUE}/vcf/{TISSUE}_samples.txt",
       samplemap = "input/{TISSUE}/map/{TISSUE}_samplemap.txt"
   log:
       err = "output/logs/{TISSUE}_samplemap.err",
       out = "output/logs/{TISSUE}_samplemap.out"
   shell:
       """
        module load R
        module load BCFtools
        bcftools query -l {input.vcf} > {output.vcfSamples}
        Rscript --vanilla scripts/makeSampleMap.R {output.vcfSamples} {output.samplemap}
        """
rule cisVar:
   input:
        expMat = "output/expression/{TISSUE}/{TISSUE}_expressionMatrix.txt.gw_phased.bed.gz",
        vcf = "input/{TISSUE}/vcf/{TISSUE}_vcf.gz",
        pair = lambda wildcards: f"ref/pairs/{wildcards.PAIRS}",
        map = "input/{TISSUE}/map/{TISSUE}_samplemap.txt"
   output:
       cisvar = "output/cisvar/{TISSUE}/{TISSUE}_{PAIRS}_pairs.txt"
   log:
       err = "output/logs/{TISSUE}_{PAIRS}_cisvar.err",
       out = "output/logs/{TISSUE}_{PAIRS}_cisvar.out"
   shell:
       """
        module load BCFtools
        python3 scripts/phaser_cis_var_ed.py --bed {input.expMat} --vcf {input.vcf} --pairs {input.pair} --map {input.map} --o {output.cisvar} --chr chr22
        """
