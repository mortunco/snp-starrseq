# Table of contents
1. [Installation](#installation)
2. [Running Example Pipeline](#running-small-example)
3. [Configuration and Outputs](#configuration-and-outputs)
4. [Visualisation](#visualisation)
5. [Publications](#publication)
6. [Contact and Support](#contant-and-support)

## Installation

### 1. Create a project directory and its sub-structures.
```
mkdir test && mkdir test/code && mkdir test/analysis && mkdir test/raw-data && mkdir test/code
cd test/code ## and enter the code directory under test
```

### 2. Create a new environment for snakemake (version >= 6.4.1) You can skip this step if you already have one.
```
conda create -c bioconda -n snakemake snakemake_env
```

### 3. Clone following repositories
```
git clone https://github.com/mortunco/snp-starrseq.git
git clone -b v0.3.4 https://github.com/vpc-ccg/calib.git
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
(calib-env)$ ./calib
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
(calib-env)$ ./consensus/calib_cons
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

## Running Small Example
### Dry run to test & Run the small test data analysis
We generated a small example dataset and corresponding configuration file (config/config.small-example.yaml) to familirize our users with our pipeline. This configuration file contains required parameters to run full pipeline in asymmetric settings. For sake of simplicity, we incorporated only single sample but defined as different samples in the configuration.

``` 
## Download small data set 
cd test/raw-data # go to the test directory that we created initially.
wget https://figshare.com/ndownloader/files/31851779 -O data.tar.gz
tar -xvf data.tar.gz
cd ..

## activate preferred snakemake env
conda activate snakemake_env 
## run the pipeline
snakemake --snakefile code/snp-starrseq/Snakemake.py -j5 --use-conda --configfile code/snp-starrseq/config/config.small-example.yaml -p -n
snakemake --snakefile code/snp-starrseq/Snakemake.py -j5 --use-conda --configfile code/snp-starrseq/config/config.small-example.yaml -p
```

## Configuration and Outputs
Our pipeline needs a config file for execution. It is recommended to create a separate config file for each analysis for reproducibility. Example config file shared in the `config/config.small-example.yaml` with descriptions. Parameters of the manuscript data is in the `config/full-calibstrict.yaml`. Mandatory/optional fields are also referred in the `config/small-example.yaml`.

Our basic file structure is as follows. There are 4 main steps;
In step 1, we cluster and generate consensus sequences from asymmetric reads.
- asym.shortlong/shortlong.cluster --> Calib cluster output from independent asymmteric sequnecing runs which contain cluster number and read information.
- asym.cons.longshort/shortlong.rX.fastq --> Calib consensus output from independent asymmteric sequnecing runs which contain consensus sequence of all clusters which =>2 members.
- asym.merged.rX.fastq --> We merge longshort/shortlongs r1 and r2 in two separate files for `2-read-matching` step.
In step 2, we match long mates of the two asymmeric runs using 24 bp barcodes. 
- (5' ~ 3bp UMI - 9 bp Fragment Sequence - ... - 9 bp Fragment Sequence - 3bp UMI ~ 3')
- asym.master-barcode-cid.txt --> Master table that contains information of each unique fragment. There are possible three outcomes, clustered means fragments are present in both (longshort/shortlong) runs. Orphan means fragment is found only in single run (one of longshort/shortlong). Problematic ones such as multiple and same same were cases that same barcodes was found by different fragments in the single run. In our benchmark this rate was <%1 therefore, we dump those reads to a file for further investigation.
In step 3, we collapsed matched long asymmetric reads and retreive enhancer fragments.
- asym/pb.collapsed-fragments.bam/fastq --> As a result of the collapsing process, this file contains sucessfully collapsed fragments reference genome alignment and reads. If pacbio method is pb prefix is used.
- asym/pb.barcode-variant-table.tsv --> This file contains all mismatches/indels on the collapsed fragments with respect to reference genome. If pacbio method is used. pb prefix is used.
- asym/pb.barcode-allele.tsv --> From the previous file, we annotate each collapsed fragment whether they support alternative allele and WT allele. Only SNP events that are located in `capture_bed` bed file were considered. If pacbio method is used. pb prefix is used.
In step 4, we quantify barcode counts from symmterical starrseq data. (optional)
- samplename[481/lib].bam --> For each symmetric sample, we initially map reads to the reference genome.
- samplename[481/lib].vis-info --> File contains enhancer fragment index, enhancer fragment position and supporting symmetrical short read name. This information is required for visualisation process.
- startposcounts.samplename[481/lib].txt --> for each file, using the alignment files we generate index (startpos:6bpUMI) for each read and quantify the occurance in the sample. These indexes should be concordant with the indexes from asymmetric section of the analysis.
In step 5, we calculate bi-allelic activity of SNPs based on the change in the enhancer fragments. 
- comparison.result-table.txt --> Containts output for each different comparison that was inputted in the config file.
In step 6, we generate two separate bedgraph files for WT and VAR alleles for those signifcant SNPs in each comparison. 
- comp_a/SNPid_chr-pos-ref-alt.WT/VAR.bedGraph --> Each SNP will get a pair of bedGraph files normalized by unique plasmid count.
- vis-done --> indicated all visualisations are succesfully completed for this comparison.

Below, the expected output directory tree for small-example with asymmetric mode.
```
$ tree small-example/
├── 1-cluster-consensus
│   ├── asym.cons.longshort.r1.fastq
│   ├── asym.cons.longshort.r2.fastq
│   ├── asym.cons.shortlong.r1.fastq
│   ├── asym.cons.shortlong.r2.fastq
│   ├── asym.longshort.cluster
│   ├── asym.merged.r1.fastq
│   ├── asym.merged.r2.fastq
│   └── asym.shortlong.cluster
├── 2-match-reads
│   ├── asym.clustered.r1.fastq
│   ├── asym.clustered.r2.fastq
│   ├── asym.master-barcode-cid.txt
│   ├── asym.orphan.fastq
│   ├── asym.problematic_multiple.interleaved.fastq
│   └── asym.problematic_samesame.interleaved.fastq
├── 3-generate-fragment-lib
│   ├── asym.barcode-allele.tsv
│   ├── asym.barcode-variant-table.tsv
│   ├── asym.collapsed-fragments.bam
│   ├── asym.collapsed-fragments.bam.bai
│   ├── asym.trimmed.clustered.r1.fastq
│   ├── asym.trimmed.clustered.r2.fastq
│   ├── pb.barcode-allele.tsv ### Pipeline was run second time to generate Pacbio based fies
│   ├── pb.barcode-variant-table.tsv ### Pipeline was run second time to generate Pacbio based fies
│   ├── pb.collapsed-fragments.bam ### Pipeline was run second time to generate Pacbio based fies
│   ├── pb.collapsed-fragments.bam.bai ### Pipeline was run second time to generate Pacbio based fies
│   └── pb.fasta
├── 4-symmetric-barcode-quantification
│   ├── 481.bam
│   ├── 481.vis-info
│   ├── 482.bam
│   ├── 482.vis-info
│   ├── 483.bam
│   ├── 483.vis-info
│   ├── lib.bam
│   ├── lib.vis-info
│   ├── startposcounts.481.txt
│   ├── startposcounts.482.txt
│   ├── startposcounts.483.txt
│   └── startposcounts.lib.txt
├── 5-bi-allelic-comparisons
│   ├── comp_a.result-table.txt
│   ├── comp_b.result-table.txt
│   ├── hg19.fa.fai -> /home/tmorova/tools/hg19.fa.fai
│   ├── hg19.genome
│   └── tmp
└── 6-vis
    └── comp_a
        ├── rs10457185_chr6-109322477-G-C.REF.bedGraph
        ├── rs10457185_chr6-109322477-G-C.VAR.bedGraph
        ├── rs2273668_chr6-109323519-G-T.REF.bedGraph
        ├── rs2273668_chr6-109323519-G-T.VAR.bedGraph
        ├── rs75675305_chr6-109324353-T-C.REF.bedGraph
        ├── rs75675305_chr6-109324353-T-C.VAR.bedGraph
        ├── rs7761290_chr6-109326621-T-G.REF.bedGraph
        ├── rs7761290_chr6-109326621-T-G.VAR.bedGraph
        └── vis-done
```

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

## Pacbio CCS mode mode
snpSTARRseq pipeline can take Pacbio CCS reads as an input. To use this mode, our pipeline expect raw Pacbio output(XXX.subreads.bam) from the sequencer to be processed by [SMRT Tools](https://www.pacb.com/support/software-downloads/) 's (release_9.0.0.92188 used in our project) CCS module with default options.

```
ccs --min-length 100 XX.subreads.bam ccs.out.bam
```

## Visualisation
6th step of our pipeline generates two separate bedGraph files for WT and VAR allales all SNPs which passes minimum enhancer fragment support threshold (Default 15). These files can then be visualised with genome browsers that are compatible with BedGraph format. Our publication's figure 2C was generated by saving 15kb window of target SNPs using the BedGraph output generated by our pipeline. 

## Publication
If you have used or influenced by our method please cite our paper --> DOI/Citation.

## Contant and Support
Feel free to open an issue on our issue section.

