# phaserPopPipe

A snakemake adaption of the phASER-POP tool by Stephane Castel.

## PhASER Expression Matrix

The first step of this pipeline and phASER-POP is to generate a sample x gene matrix with all gene-level haplotypes. 

### Requirements
* `gene_ae.txt` files that are created with phASER

* Gencode annotation `bed file`

##

1.  Clone this directory into your working directory:
    ``` bash
    git clone https://github.com/mbond0718/phaserPopPipe.git .
    ```

2. Edit the tab-separated samplesheet_expr.txt file with the names of your tissue and sample IDs. This pipeline was written specifically to process dGTEx data, so it groups samples into the tissues that they belong to. Additional metadata can be included here. 
    | Tissue  | SampleID  | FullID   |
    |---------|-----------|----------|
    | LUNG    | D001      | D001-123 |
    | LUNG    | D002      | D002-123 | 
    | MUSCLE  | D001      | D001-456 |
    | MUSCLE  | D002      | D002-456 |
    | COLON   | D001      | D001-789 |
    | COLON   | D002      | D002-789 |

3. Edit `config/config.yaml` file for your system/experiment. 
    ```yaml
    ## Path to sample sheet for this workflow and the qtl (see below)
    samplesheet: 'samplesheet_expr.txt'
    samplesheet_qtl: 'samplesheet_cisvar.txt'

    ## Other parameters
    ref: "/gpfs/commons/groups/lappalainen_lab/mbond/work/phaserPop/ref/gencode.v25.annotation.bed"
    aseDir: "/gpfs/commons/groups/lappalainen_lab/mbond/work/asePipe_V3/output/ase/"
    outpath: "/gpfs/commons/groups/lappalainen_lab/mbond/work/phaserPopPipe/"
    ```

4.  Submit to `SLURM` with sbatch script:
    ```bash
    sbatch launchExpression.sh
    ```
 `launchExpression.sh` will submit a long-running, low-resource job script that will spawn other jobs in the pipeline as dependencies are fulfilled. For systems other than slurm, edit `launchExpression.sh` to create a submission script for your HPC's job scheduler.

After running these steps the pipeline will produce the following files:
- `output/expression/{TISSUE}/{TISSUE}_expressionMatrix.txt.bed.gz`
- `expression/{TISSUE}/{TISSUE}_expressionMatrix.txt.bed.gz.tbi`
- `output/expression/{TISSUE}/{TISSUE}_expressionMatrix.txt.gw_phased.bed.gz`
- `output/expression/{TISSUE}/{TISSUE}_expressionMatrix.txt.gw_phased.bed.gz.tbi`


## phASER Cis Var

This second part of the phASER-POP pipeline that estimates effect sizes for a list of variants using aggregated haplotype level data from phASER and the genotypes of your samples. This pipeline was developed specifically to process dGTEx data, so it uses adult GTEx eQTLs as the variants and compares all adult GTEx tissues to dGTEx tissues in parallel. 

### Requirements
* expressionMatrix.txt.gw_phased.bed.gz (previous workflow)
* `vcf` with sample genotypes
* `variant list` to calculate effect sizes (for this workflow specifically, I used the GTEx v11 eQTL pairs `GTEx_Analysis_v11_eQTL.tar` from the GTEx portal)

##

1. Edit the tab-separated samplesheet_cisvar.txt file with the names of your tissues and pairs files to use. Since this workflow was written for dGTEx, it uses all of the samples for an given tissue, but this can be edited to be run on an individual level.
    | Tissue  | Pairs     | 
    |---------|-----------|
    | LUNG    | lung_pairs.txt      | 
    | LUNG    | blood_pairs.txt      |
    | MUSCLE  | muscle_pairs.txt |
    | MUSCLE  | blood_pairs.txt      |
    | COLON   | colon_pairs.txt      |
    | COLON   | blood_pairs.txt      |

2.  Submit to `SLURM` with sbatch script:
    ```bash
    sbatch launchCisVar.sh
    ```
 `launchCisVar.sh` will submit a long-running, low-resource job script that will spawn other jobs in the pipeline as dependencies are fulfilled. For systems other than slurm, edit `launchExpression.sh` to create a submission script for your HPC's job scheduler.

After running these steps the pipeline will produce the following files:
- `input/{TISSUE}/vcf/{TISSUE}_vcf.gz`
- `input/{TISSUE}/vcf/{TISSUE}_vcf.gz.tbi`
- `input/{TISSUE}/vcf/{TISSUE}_samples.txt`
- `input/{TISSUE}/map/{TISSUE}_samplemap.txt`
- `output/cisvar/{TISSUE}/{TISSUE}_{QTL}_pairs.txt`

Output directory structure:
```
slurm-{jobid}.out
working_directory/
├── input/
│   └── {group}/
│       ├── map/
│       │   └── samples.txt
│       └── vcf/
│           ├── {group}.vcf.gz
│           └── {group}.vcf.gz.tbi
└── output/
    ├── expression/
    │   └── {group}/
    │       └── expressionMatrix.txt.gw_phased.bed.gz
    └── cisvar/
        └── {group}/
            └── {group}_{pairs}.txt
```
### References
- [phASER-POP](https://github.com/secastel/phaser/tree/master/phaser_pop)
