
### snakemake checkes config parameters to be defined despite they are used in the run. 
### Following lines check and fill those variable if they are not used to prevent error.

if list(config["fragment_retreival"].keys())[0] == "pb": ### means this is a pacbio based run
  if config["smrttools_full_path"] == "" or "smrttools_full_path" not in config:
    raise ValueError('Reconstruction of fragments are based on Pacbio reads but SMRTtools path (smrttools_full_path) is not defined')
  if config["parts"] == "" or "parts" not in config:
    raise ValueError('Reconstruction of fragments are based on Pacbio reads but number of parallel parts (parts) is not defined')
  if config["longread_adaptor"] == "" or "longread_adaptor" not in config:
    raise ValueError('Reconstruction of fragments are based on Pacbio reads but long reads adaptors (longread_adaptor) is not defined')
else: ### means this is an asymmetrical read based run.
  config["smrttools_full_path"]=""
  config["parts"]=10
  config["longread_adaptor"]=""
  if config["calib_params"] == "" or "calib_params" not in config: 
    raise ValueError('Reconstruction of fragments are based on asymmetrical reads but calib (calib_params) clustering parameters are not defined')
    if config["longread_validation"] == "" or "longread_validation" not in config:
      config["longread_validation"]== ""
      print(f'Reconstructed fragment confirmation with respect to Pacbio CCS wont be executed.')
    
def get_replicates(wildcards):
  if list(config["fragment_retreival"].keys())[0] == "pb":
    return(config["fragment_retreival"]["pb"]["reads"])
  else:
    return(expand("{x}",x=config["fragment_retreival"]["asym"][wildcards.sample_name][wildcards.direction]))

def get_snp_annotation(wildcards):
  if config["snp_annotation"] == "" or "snp_annotation" not in config:
    return("noannot")
  else:
    return(config["snp_annotation"])

def get_min_count(wildcards):
  if config["bi_allelic_NBR_mincount"] == "" or "bi_allelic_NBR_mincount" not in config:
    return(15)
  else:
    return(config["bi_allelic_NBR_mincount"])

def get_stat_biallele_barcode_allele(wildcards):
  return('analysis/{run_name}/3-generate-fragment-lib/{seq_type}.barcode-allele.tsv'.format(run_name=config["run_name"],seq_type=list(config["fragment_retreival"].keys())[0]) )

def get_stat_biallele_mrna(wildcards):
  return(expand("analysis/{rn}/4-symmetric-barcode-quantification/startposcounts.{x}.txt",rn=config["run_name"],x=config["bi_allelic_comparisons"][wildcards.comp_name]["mrna"]))

def get_stat_biallele_mrna2(wildcards): ## this function is here to make snakemake work but output with better formatting for argparse 
  return(",".join(expand("analysis/{rn}/4-symmetric-barcode-quantification/startposcounts.{x}.txt",rn=config["run_name"],x=config["bi_allelic_comparisons"][wildcards.comp_name]["mrna"])))

def get_stat_biallele_mrna_vis(wildcards):
  return(expand("analysis/{rn}/4-symmetric-barcode-quantification/{x}.vis-info",rn=config["run_name"],x=config["bi_allelic_comparisons"][wildcards.comp_name]["mrna"]))

def get_stat_biallele_mrna_vis2(wildcards): ## this function is here to trick snakemake work but output with better formatting for argparse.
  return(",".join(expand("analysis/{rn}/4-symmetric-barcode-quantification/{x}.vis-info",rn=config["run_name"],x=config["bi_allelic_comparisons"][wildcards.comp_name]["mrna"])))

def get_stat_biallele_dna(wildcards):
  return(expand("analysis/{rn}/4-symmetric-barcode-quantification/startposcounts.{x}.txt",rn=config["run_name"],x=config["bi_allelic_comparisons"][wildcards.comp_name]["dna"]))

def get_mrna_r1(wildcards):
  return(expand("{x}",x=config["quantification"][wildcards.mrna_sample]["r1"]))
def get_mrna_r2(wildcards):
  return(expand("{x}",x=config["quantification"][wildcards.mrna_sample]["r2"]))

def generate_results():
  results=[]
  run_name=config["run_name"]
  if "asym" in config["fragment_retreival"]:
    results.append(f'analysis/{run_name}/3-generate-fragment-lib/asym.barcode-allele.tsv')
    results.append(f'analysis/{run_name}/3-generate-fragment-lib/asym.collapsed-fragments.bam')
  else:
    results.append(f'analysis/{run_name}/3-generate-fragment-lib/pb.barcode-allele.tsv')
    results.append(f'analysis/{run_name}/3-generate-fragment-lib/pb.collapsed-fragments.bam')
  for i in list(config["quantification"].keys()):
      results.append(f'analysis/{run_name}/4-symmetric-barcode-quantification/startposcounts.{i}.txt')
  if config["longread_validation"] != "" and "longread_validationongread_samples" in config:
    results.append(f'analysis/{run_name}/misc-longread/longread-counts.txt')
  if config["bi_allelic_comparisons"] != "" and "bi_allelic_comparisons" in config:
    for comp_name,values in config["bi_allelic_comparisons"].items():
      results.append(f'analysis/{run_name}/5-bi-allelic-comparisons/{comp_name}.result-table.txt')
      results.append(f'analysis/{run_name}/6-vis/{comp_name}/vis-done')
  return results

rule all:
  input: generate_results()

rule asym_prep_merge_replicates:
  input:
    get_replicates
  output:
    temp("analysis/{run_name}/1-cluster-consensus/asym.{sample_name}.{direction}.fastq")
  conda:
    "env/env.seqtk.yaml"
  threads:1
  shell:
    """
    zcat {input} > {output}
    """
rule asym_trim_for_calib:
  input:
    "analysis/{run_name}/1-cluster-consensus/asym.{sample_name}.{direction}.fastq"
  output:
    temp("analysis/{run_name}/1-cluster-consensus/asym.trimmed.{sample_name}.{direction}.fastq")
  conda:
    "env/env.seqtk.yaml"
  threads:1
  shell:
    """
    seqtk trimfq -L 75 {input} > {output}
    """

rule asym_calib_cluster:
  input:
    r1="analysis/{run_name}/1-cluster-consensus/asym.trimmed.{sample_name}.r1.fastq",
    r2="analysis/{run_name}/1-cluster-consensus/asym.trimmed.{sample_name}.r2.fastq"
  output:
    "analysis/{run_name}/1-cluster-consensus/asym.{sample_name}.cluster"
  params:
    config["calib_params"]
  log: 
    "logs/{run_name}/1-calib-cluster.asym.{sample_name}.log"
  threads:4
  shell:
    """
    ./code/calib/calib/calib -f {input.r1} -r {input.r2} -o analysis/{wildcards.run_name}/1-cluster-consensus/asym.{wildcards.sample_name}. {params} 2> {log}
    """

rule asym_calib_cons:
  input:
    r1="analysis/{run_name}/1-cluster-consensus/asym.{sample_name}.r1.fastq",
    r2="analysis/{run_name}/1-cluster-consensus/asym.{sample_name}.r2.fastq",
    cluster="analysis/{run_name}/1-cluster-consensus/asym.{sample_name}.cluster"
  output:
    r1_cons_tmp=temp("analysis/{run_name}/1-cluster-consensus/asym.constemp.{sample_name}.r1.fastq"),
    r2_cons_tmp=temp("analysis/{run_name}/1-cluster-consensus/asym.constemp.{sample_name}.r2.fastq"),
    r1_cons_msa_tmp=temp("analysis/{run_name}/1-cluster-consensus/asym.constemp.{sample_name}.r1.msa"),
    r2_cons_msa_tmp=temp("analysis/{run_name}/1-cluster-consensus/asym.constemp.{sample_name}.r2.msa")
  log: 
    "logs/{run_name}/2-calib-cons.asym.{sample_name}.log"
  threads: 4
  shell:
    """
    ./code/calib/calib/consensus/calib_cons -q {input.r1} -q {input.r2} -o analysis/{wildcards.run_name}/1-cluster-consensus/asym.constemp.{wildcards.sample_name}.r1 analysis/{wildcards.run_name}/1-cluster-consensus/asym.constemp.{wildcards.sample_name}.r2 -c {input.cluster} 2> {log}
    """

rule asym_label_reads:
  input:
    r1="analysis/{run_name}/1-cluster-consensus/asym.constemp.{sample_name}.{direction}.fastq",
  output:
    r1="analysis/{run_name}/1-cluster-consensus/asym.cons.{sample_name}.{direction}.fastq",
  threads:1
  shell:
    """
    cat {input.r1}  | awk 'BEGIN {{OFS="\t"}} {{if(NR % 4 == 1) {{print $1_\"{wildcards.sample_name}\",$2;}} else print $0 }};' > {output.r1}
    """


rule asym_gather_reads:
  input:
    a="analysis/{run_name}/1-cluster-consensus/asym.cons.longshort.{direction}.fastq",
    b="analysis/{run_name}/1-cluster-consensus/asym.cons.shortlong.{direction}.fastq"
  output:
    r1="analysis/{run_name}/1-cluster-consensus/asym.merged.{direction}.fastq"
  threads:1
  shell:
    """
    cat {input.a} {input.b}  > {output.r1}
    """

rule asym_match_reads:
  input:
    r1="analysis/{run_name}/1-cluster-consensus/asym.merged.r1.fastq",
    r2="analysis/{run_name}/1-cluster-consensus/asym.merged.r2.fastq"
  output:
    orphan="analysis/{run_name}/2-match-reads/asym.orphan.fastq",
    clustered_r1="analysis/{run_name}/2-match-reads/asym.clustered.r1.fastq",
    clustered_r2="analysis/{run_name}/2-match-reads/asym.clustered.r2.fastq",
    problematic_samesame="analysis/{run_name}/2-match-reads/asym.problematic_samesame.interleaved.fastq",
    problematic_multiple="analysis/{run_name}/2-match-reads/asym.problematic_multiple.interleaved.fastq",
    master_file="analysis/{run_name}/2-match-reads/asym.master-barcode-cid.txt"
  threads:1
  shell:
    """
    python code/snp-starrseq/read_id_matcher.py --f {input.r1} --r {input.r2} --p analysis/{wildcards.run_name}/2-match-reads/asym.
    """

rule asym_pre_collapse_trim:
  input:
    clust_r1="analysis/{run_name}/2-match-reads/asym.clustered.r1.fastq",
    clust_r2="analysis/{run_name}/2-match-reads/asym.clustered.r2.fastq"
  output:
    trim_r1=temp("analysis/{run_name}/3-generate-fragment-lib/asym.trimmed.clustered.r1.fastq"),
    trim_r2=temp("analysis/{run_name}/3-generate-fragment-lib/asym.trimmed.clustered.r2.fastq"),
  conda:
    "env/env.seqtk.yaml"
  threads:1
  shell:
    """
      seqtk trimfq -q 0.05 {input.clust_r1} > {output.trim_r1}
      seqtk trimfq -q 0.05 {input.clust_r2} > {output.trim_r2}
    """

rule asym_collapse_fragments:
  input:
    trim_r1=temp("analysis/{run_name}/3-generate-fragment-lib/asym.trimmed.clustered.r1.fastq"),
    trim_r2=temp("analysis/{run_name}/3-generate-fragment-lib/asym.trimmed.clustered.r2.fastq"),
  output:
    collapsed=temp("analysis/{run_name}/3-generate-fragment-lib/asym.collapsed-fragments.fastq")
  conda: 
    "env/env.collapse.yaml"
  log: 
    "logs/{run_name}/3-bbmerge.asym.log"
  threads:1
  shell:
    """
    bbmerge.sh -Xmx100G ecco=t merge=t in1={input.trim_r1} in2={input.trim_r2} out={output.collapsed} 2>{log}
    """

rule asym_map_collapsed_to_ref:
  input: 
    "analysis/{run_name}/3-generate-fragment-lib/asym.collapsed-fragments.fastq"
  output:
    "analysis/{run_name}/3-generate-fragment-lib/asym.collapsed-fragments.bam"
  conda:
    "env/env.mapping.yaml"
  log: 
    "logs/{run_name}/3-map-fragments.asym.log"
  params:
    config["bwa_ref"]
  threads:16
  shell:
    """
    bwa mem -t{threads} {params} {input} | samtools view -Sb | samtools sort -@{threads} -m10G> {output} && samtools index {output} 2> {log}
    """

rule pacbio_longread_css:
  input:
    get_replicates
  output:
    "analysis/{run_name}/3-generate-fragment-lib/pb.bam"
  conda:
    "env/env.pacbio.yaml"
  params:
    config["smrttools_full_path"]
  shell:
    """
    {params} --min-length 100 {input} {output}
    """
rule pacbio_extract_reads:
  input:
    "analysis/{run_name}/3-generate-fragment-lib/pb.bam"
  output:
    temp("analysis/{run_name}/3-generate-fragment-lib/pb.fasta")
  conda:
    "env/env.pacbio.yaml"
  shell:
    """
    samtools fastq -@10 -n {input} |  awk 'NR % 4 == 1 || NR % 4 == 2' | sed 's/@/>/g' -  > {output}
    """
rule pacbio_fasta_split:
  input:
    "analysis/{run_name}/3-generate-fragment-lib/pb.fasta"
  output:
    temp(expand("analysis/{{run_name}}/3-generate-fragment-lib/pb.fasta{cnt}",cnt=range(config['parts'])))
  conda:
    "env/env.pacbio.yaml"
  shell:
    """
    split {input} pb.fasta -d -a1 -n10
    """
rule pacbio_mrsfast_index:
  input:
    "analysis/{run_name}/3-generate-fragment-lib/pb.fasta{cnt}"
  output:
    temp("analysis/{run_name}/3-generate-fragment-lib/pb.sam{cnt}")
  conda:
    "env/env.pacbio.yaml"
  params:
    config["longread_adaptor"]
  shell:
    """
    mrsfast --index {input}
    mrsfast --search {input} --seq {params} -e2 -o {output}
    """
rule pacbio_merge_extract:
  input:
    split=expand("analysis/{{run_name}}/3-generate-fragment-lib/pb.sam{cnt}",cnt=range(config['parts'])),
    full="analysis/{run_name}/3-generate-fragment-lib/pb.fasta"
  output:
    merged=temp("analysis/{run_name}/3-generate-fragment-lib/pb.merged.bam"),
    extract="analysis/{run_name}/3-generate-fragment-lib/withoutadapter.fasta"
  conda:
    "env/env.pacbio.yaml"
  params:
      param= "-m8G"
  threads:
      16
  shell:
    """
    samtools merge {output} {input.split} -@ {threads} {params.param}
    samtools view {output} | python ../sep-read-adapters.py -f {input.full} -s - -a 33 -o analysis/{wildcards.run_name}/3-generate-fragment-lib/pb
    """
rule pacbio_map_markdup:
  input:
    "analysis/{run_name}/3-generate-fragment-lib/withoutadapter.fasta"
  output:
    temp=temp("analysis/{run_name}/3-generate-fragment-lib/temp.withoutadapter.fasta"),
    temp0=temp("analysis/{run_name}/3-generate-fragment-lib/temp.withoutadapter.all0dir.sam"),
    temp0bam=temp("analysis/{run_name}/3-generate-fragment-lib/temp.withoutadapter.all0dir.bam"),
    final="analysis/{run_name}/3-generate-fragment-lib/pb.collapsed-fragments.bam",
  conda:
    "env/env.pacbio.yaml"
  shell:
    """
    minimap2 -t20 -ax asm20 -c ~/tools/hg19.fa {input} | samtools view -Sb | samtools sort -@10 -m10G > {output.temp};
    samtools view -F 2052 {output.temp} | awk 'BEGIN {{OFS="\t"}} {{print $1,0,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}}' > {output.temp0};
    cat <(samtools view -H {output.temp0}) {output.temp0} | samtools view -Sb - | samtools sort -@10 -m10G > {output.temp0bam};
    picard MarkDuplicates I={output.temp0bam} O={output.final} M=pacbio.markdup.metrics.txt VALIDATION_STRINGENCY=LENIENT;
    samtools index {output.final}
    """

rule create_variant_table:
  input:
    "analysis/{run_name}/3-generate-fragment-lib/{seq_type}.collapsed-fragments.bam"
  output:
    "analysis/{run_name}/3-generate-fragment-lib/{seq_type}.barcode-variant-table.tsv"
  conda:
    "env/env.variant-db.yaml"
  threads:1
  shell:
    """
    python code/snp-starrseq/create-barcode-variant-db.py --bam_input {input} > {output}
    """

rule create_mutation_database:
  input:
    bam="analysis/{run_name}/3-generate-fragment-lib/{seq_type}.collapsed-fragments.bam",
    variant_table="analysis/{run_name}/3-generate-fragment-lib/{seq_type}.barcode-variant-table.tsv"
  output:
    "analysis/{run_name}/3-generate-fragment-lib/{seq_type}.barcode-allele.tsv"
  conda:
    "env/env.variant-db.yaml"
  params:
    bed=config["capture_bed"]
  log: 
    "logs/{run_name}/{seq_type}.3-mutation-database.log"
  threads:1
  shell:
    """
    python code/snp-starrseq/dev-segment-variant-matrix.py --bam_input {input.bam} --var_db {input.variant_table} --csv_out {output} --mat_out False --bed_file {params.bed} 2>{log}
    """
rule mrna_quantification_startpos:
  input:
    r1=get_mrna_r1,
    r2=get_mrna_r2,
  output:
    countx="analysis/{run_name}/4-symmetric-barcode-quantification/startposcounts.{mrna_sample}.txt",
    bam="analysis/{run_name}/4-symmetric-barcode-quantification/{mrna_sample}.bam",
    vis="analysis/{run_name}/4-symmetric-barcode-quantification/{mrna_sample}.vis-info"
  conda:
    "env/env.mapping.yaml"	
  params:
    config["bwa_ref"]
  threads:16
  shell:
    """
    bwa mem -t{threads} {params} {input.r1} {input.r2} | samtools view -Sb > {output.bam};
    samtools view -f2 -F2048 {output.bam} | python code/snp-starrseq/create-umi-directory-paired.py - > {output.vis};
    cut -f 1 {output.vis} | sort -T . --parallel={threads} | uniq -c | awk '{{print $1,$2}}' > {output.countx}
    """

rule stat_comparison:
  input:
    barcode_allele=get_stat_biallele_barcode_allele, 
    dna=get_stat_biallele_dna,
    mrna=get_stat_biallele_mrna
  output:
    "analysis/{run_name}/5-bi-allelic-comparisons/{comp_name}.result-table.txt"
  params:
    annotation=get_snp_annotation,
    min_count=get_min_count,
    mrna_formatted=get_stat_biallele_mrna2
  conda:
    "env/env.NBR.yaml"
  shell:
    """
    Rscript ./code/snp-starrseq/bi-allelic-NBR.R \
      --dna={input.dna} \
      --mrna={params.mrna_formatted} \
      --barcode_allele={input.barcode_allele} \
      --output={output} \
      --annotation={params.annotation} \
      --min_count={params.min_count}
    """

rule vis_data:
  input:
    barcode_allele=get_stat_biallele_barcode_allele,
    mrna=get_stat_biallele_mrna_vis,
    restable="analysis/{run_name}/5-bi-allelic-comparisons/{comp_name}.result-table.txt"
  output:
    "analysis/{run_name}/6-vis/{comp_name}/vis-done"
  params:
    genome=config["bwa_ref"],
    mrna_formatted=get_stat_biallele_mrna_vis2
  conda:
    "env/env.vis.yaml"
  shell:
    """
    python code/snp-starrseq/generate-tracks-bi-allelic.py \
      -i {params.mrna_formatted} \
      -r {input.restable} \
      -b {input.barcode_allele} \
      -g {params.genome} -o analysis/{wildcards.run_name}/6-vis/{wildcards.comp_name} > {output}  
    """

rule misc_longread_quantification_startpos:
  input:
    config["longread_validation"]
  output:
    "analysis/{run_name}/misc/longread-counts.txt"
  conda:
    "env/env.mapping.yaml"	
  threads: 16
  shell:
    """
    samtools view -q2 -F2052 {input} | python code/snp-starrseq/create-umi-directory-longread.py | cut -f 1 | sort -T . --parallel={threads} | uniq -c  | awk '{{print $1,$2}}' > {output}
    """
