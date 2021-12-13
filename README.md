# Table of contents
1. [Installation](#installation)
2. [Running Example Pipeline](#running-pamir-example)
3. [Configuration](#configuration)
4. [Output File Definitions](#output-file-definitions)
5. [Visualisation](#visualisation)
6. [Publications](#publication)
7. [Contant and Support](#contant-and-support)

## Installation

### 1. Create a project directory and its sub-structures.
```
mkdir test && mkdir test/code && mkdir test/analysis && mkdir test/raw-data && mkdir test/code
cd test/code ## and enter the code directory under test
```

### 2. Clone following repositories
```
git clone https://github.com/mortunco/snp-starrseq.git
git clone -b v0.3.4 https://github.com/vpc-ccg/calib.git
```

### 3. Copy example data to raw-data and untar 
```
cp data/small-realdata.tar.gz ../../raw-data/
tar xvf small-realdata.tar.gz
```

### 4. Create environment for calib
```
conda env create --file snp-starrseq/env/env.calib.yaml
exit # to reset variables.
```

### 5. Build Calib from source. 
Our pipeline requires Calib for asymmetric sequencing fragmant reconstruction. 
```
exit ### Start a new terminal 
cd go/to/your/test
conda activate calib-v3.4
export CPATH=${CONDA_PREFIX}/include
make ## compiles calib
make -C consensus/ # compiles calib's consensus module
```

Once calib is compiled we dont need this environment. Calib and consensus should return the following.
```
(calib-x) [tmorova@linuxsrv006 calib]$ ./calib
Combined barcode lengths must be a positive integer and each mate barcode length must be non-negative! Note if both mates have the same barcode length you can use -l/--barcode-length parameter instead.
Calib: Clustering without alignment using LSH and MinHashing of barcoded reads
Usage: calib [--PARAMETER VALUE]
Example: calib -f R1.fastq -r R2.fastq -o my_out. -e 1 -l 8 -m 5 -t 2 -k 4 --silent
Calib's paramters arguments:
        -f    --input-forward                   (type: string;   REQUIRED paramter)
        -r    --input-reverse                   (type: string;   REQUIRED paramter)
        -o    --output-prefix                   (type: string;   REQUIRED paramter)
        -s    --silent                          (type: no value; default: unset)
        -q    --no-sort                         (type: no value; default:  unset)
        -g    --gzip-input                      (type: no value; default:  unset)
        -l    --barcode-length                  (type: int;      REQUIRED paramter unless -l1 and -l2 are provided)
        -l1   --barcode-length-1                (type: int;      REQUIRED paramter unless -l is provided)
        -l2   --barcode-length-2                (type: int;      REQUIRED paramter unless -l is provided)
        -p    --ignored-sequence-prefix-length  (type: int;      default: 0)
        -m    --minimizer-count                 (type: int;      default: Depends on observed read length;)
        -k    --kmer-size                       (type: int;      default: Depends on observed read length;)
        -e    --error-tolerance                 (type: int;      default: Depends on observed read length;)
        -t    --minimizer-threshold             (type: int;      default: Depends on observed read length;)
        -c    --threads                         (type: int;      default: 1)
        -h    --help
```

```
(calib-x) [tmorova@linuxsrv006 calib]$ ./consensus/calib_cons
No cluster filename was passed.
Calib Consensus: Generating consensus sequence from Calib clusters.
Usage: calib_cons [--PARAMETER VALUE]
Example 1: calib_cons -t 8 -c input.cluster -q 1.fastq 2.fastq -o 1.out 2.out
Example 2: calib_cons -q 1.fastq -q 2.fastq -o 1.out 2.out -c input.cluster
Calib's paramters arguments:
  -q  --fastq                    (type: space separated string list;
                                    REQUIRED paramter;
                                    can be set multiple times like in Example 2)
  -o  --output-prefix            (type: space separated string list;
                                    REQUIRED paramter;
                                    can be set multiple times like in Example 2;
                                    must be same size as fastq list)
  -c  --cluster                  (string;
                                    REQUIRED paramter)
  -t  --threads                  (positive integer;
                                    default: 4)
  -m  --min-reads-per-cluster    (positive integer;
                                    default: 2)
  -h  --help
```

### Create a new environment for snakemake. You can skip this step if you have one.
```
conda create -c bioconda -n snakemake snakemake_env
```

## Running Small Example
### Dry run to test & Run the small test data analysis 
```
conda activate snakemake_env 
cd go/to/your/test
snakemake --snakefile code/snp-starrseq/Snakemake.py -j5 --use-conda --configfile code/snp-starrseq/config/config.small-example.yaml -p -n
snakemake --snakefile code/snp-starrseq/Snakemake.py -j5 --use-conda --configfile code/snp-starrseq/config/config.small-example.yaml -p

```
## Configuration
Our pipeline needs a config file for execution. It is recommended to create a separate config file for each analysis for reproducibility. Following options were are in configuration file. Argument names with * are not optional, they are mandatory. Please see the config.small-example.yaml for example run. Manuscript data could be reproduced by config.full-calibstrict.yaml

| Name               | Description                                                   |
|--------------------|---------------------------------------------------------------|
| run_name*           | name of the run                                               |
| bwa_ref*            | Path to BWA index files                                       |
| calib_params*       | Calib's clustering module parameters                          |
| fragment_retreival*       | Path to asymmetric samples.                                   |
|                    | Fastq files must be in zipped format.                         |
|                    | In case of replicates, use "," to merge.                      |
|                    | For example,                                                  |
|                    | r1: ["a.r1.fastq","b.r2.fastq"]                               |
|                    | r2: ["a.r2.fastq","b.r2.fastq"]                               |
| capture_bed*        | STARRseq Capture Regions in BED format.                       |
| symmetric _samples | Path to symmteric samples.                                    |
|                    | Exteded options are same as asymmetric.                       |
| longread_samples   | Path to longread bam file.                                    |
|                    | Must be mapping to same the reference genome as other samples |

## Output File Definitions
Our basic file structure is as follows. There are 4 main steps;
- In step 1, we cluster and generate consensus sequences from asymmetric reads.
- - longshort/shortlong.cluster --> Calib cluster output which contains cluster number and read information. 
- - conx.longshort/shortlong.rX.fastq --> Calib consensus output which contains consensus sequence of all clusters which =>2 members.
- - merged.rX.fastq --> We merge longshort/shortlongs r1 and r2 in two separate files for `2-read-matching` step.
- In step 2, we match long mates of the two asymmeric runs using 24 bp barcodes. 
- - (5' ~ 3bp UMI - 9 bp Fragment Sequence - ... - 9 bp Fragment Sequence - 3bp UMI ~ 3')
- - master-barcode-cid.txt --> Master table that contains information of each unique fragment. There are possible three outcomes, clustered means fragments are present in both (longshort/shortlong) runs. Orphan means fragment is found only in single run (one of longshort/shortlong). Problematic ones such as multiple and same same were cases that same barcodes was found by different fragments in the single run. In our benchmark this rate was <%1 therefore, we dump those reads to a file for further investigation.
- In step 3, we collapsed matched long asymmetric reads and retreive enhancer fragments.
- - trimmed.clustered.rX.fastq --> Bad quality sequences and N bases were trimmed from clustered reads (check previous section).
- - collapsed-fragments.bam/fastq --> As a result of the collapsing process, this file contains sucessfully collapsed fragments reference genome alignment and reads. 
- - barcode-variant-table.tsv --> This file contains all mismatches/indels on the collapsed fragments with respect to reference genome.
- - barcode-allele.tsv --> From the previous file, we annotate each collapsed fragment whether they support alternative allele and WT allele. Only SNP events that are located in `capture_bed` bed file were considered.
- In step 4, we quantify barcode counts from symmterical starrseq data. (optional)
- - samplename[481/lib].bam --> For each symmetric sample, we initially map reads to the reference genome.
- - samplename[481/lib].sorted.bam --> sorted by position.
- - startposcounts.samplename[481/lib].txt --> for each file, using the alignment files we generate index (startpos:6bpUMI) for each read and quantify the occurance in the sample. These indexes should be concordant with the indexes from asymmetric section of the analysis.
- In step 5, unique fragments library can also be retreived from long-read sequences (Pacbio?).
- - longread-counts.txt --> For validation purposes, we sequenced our plasmid library with Pacbio CCS reads. Similar to symmetrical samples, we generate index for each alignment (which > mapq2 and not supplementary) and quantify the number of the occurances.

```
$ tree analysis/
analysis/
└── small-example
    ├── 1-cluster-consensus
    │   ├── cons.longshort.r1.fastq
    │   ├── cons.longshort.r2.fastq
    │   ├── cons.shortlong.r1.fastq
    │   ├── cons.shortlong.r2.fastq
    │   ├── longshort.cluster
    │   ├── merged.r1.fastq
    │   ├── merged.r2.fastq
    │   └── shortlong.cluster
    ├── 2-match-reads
    │   ├── clustered.r1.fastq
    │   ├── clustered.r2.fastq
    │   ├── master-barcode-cid.txt
    │   ├── orphan.fastq
    │   ├── problematic_multiple.interleaved.fastq
    │   └── problematic_samesame.interleaved.fastq
    ├── 3-generate-fragment-lib
    │   ├── barcode-allele.tsv
    │   ├── barcode-variant-table.tsv
    │   ├── collapsed-fragments.bam
    │   ├── collapsed-fragments.bam.bai
    │   ├── trimmed.clustered.r1.fastq
    │   └── trimmed.clustered.r2.fastq
    ├── 4-symmetric-barcode-quantification
    │   ├── 481.bam
    │   ├── 481.sorted.bam
    │   ├── 481.sorted.bam.bai
    │   ├── lib.bam
    │   ├── lib.sorted.bam
    │   ├── lib.sorted.bam.bai
    │   ├── startposcounts.481.txt
    │   └── startposcounts.lib.txt
    └── 5-longread
        └── longread-counts.txt
```
## Visualisation
Lorem ipsum

## Annotation file for RSids.
5th step of our pipeline calculates bi-allelic activity of each SNP with our novel NBR method. To include RSid (rsXXX) for the events our analysis code requires annotation file which is consisted two column such as 1-base genomic position with chr:pos (**no chr**) and RS snp id (rsXXX). We shared dbsnp150 common VCF annotation for hg19/grch37 and grch38 genomes in our repository. But any user specific VCF could be used for this procedure if other annotation is neccasary.

```
# for example, common VCF file from dbsnp150
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/00-common_all.vcf.gz -O common_all.dbsnp150.grch38.vcf.gz
zcat common_all.dbsnp150.grch38.vcf.gz | grep -v "#" | cut -f 1-3 | awk '{print $1 ":" $2 "\t" $3}' > dbsnp150.grch38.snp

$head dbsnp150.grch38.snp
1:10177 rs367896724
1:10352 rs555500075
1:10352 rs145072688
1:10616 rs376342519
1:10642 rs558604819
1:11008 rs575272151
1:11012 rs544419019
1:11063 rs561109771
1:13110 rs540538026
```

## Publication
TBD
## Contant and Support
Lorem Ipsum

