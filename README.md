# Table of contents
1. [Installation](#Installation)
2. [Running the Pipeline](#running-pamir)
2. [Outputs](#Outputs)

## Installation

### create a project directory and its sub-structures.
```
mkdir test && mkdir test/code && mkdir test/analysis && mkdir test/raw-data && mkdir test/code
cd test/code ## and enter the code directory under test
```

### clone neccesary repositories
```
git clone https://github.com/mortunco/snp-starrseq.git
git clone -b v0.3.4 https://github.com/vpc-ccg/calib.git
```

### copy example data to raw-data and untar 
```
cp data/small-realdata.tar.gz ../../raw-data/
tar xvf small-realdata.tar.gz
```

### create environment for calib
```
conda env create --file snp-starrseq/env/env.calib.yaml
exit # to reset variables.
```
### start new terminal ###
```
cd go/to/your/test`
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

## Running Pamir
### Dry run to test & Run the small test data analysis 
```
conda activate snakemake_env 
cd go/to/your/test
snakemake --snakefile code/snp-starrseq/Snakemake.py -j5 --use-conda --configfile code/snp-starrseq/config/config.small-example.yaml -p -n
snakemake --snakefile code/snp-starrseq/Snakemake.py -j5 --use-conda --configfile code/snp-starrseq/config/config.small-example.yaml -p
```

## Outputs
Our basic file structure is as follows. There are 4 main steps;
- In step 1, we cluster and generate consensus sequences from asymmetric reads.
- In step 2, we match long mates of the two asymmeric runs using 24 bp barcodes. 
- - (5' ~ 3bp UMI - 9 bp Fragment Sequence - ... - 9 bp Fragment Sequence - 3bp UMI ~ 3')
- In step 3, we collapsed matched long asymmetric reads and retreive enhancer fragments.
- In step 4, we quantify barcode counts from symmterical starrseq data. (optional)
- In step 5, unique fragments library can also be retreived from long-read sequences (Pacbio?).

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

