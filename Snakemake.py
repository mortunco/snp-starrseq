def get_replicates(wildcards):
	#print(expand("raw-data/{x}",x=config["samples"][wildcards.sample_name][wildcards.direction] ))
	return(expand("{x}",x=config["dna_samples"][wildcards.sample_name][wildcards.direction] ))
def get_mrna_r1(wildcards):
	#print(expand("raw-data/{x}",x=config["samples"][wildcards.sample_name][wildcards.direction] ))
	return(expand("{x}",x=config["mrna_samples"][wildcards.mrna_sample]["r1"] ))
def get_mrna_r2(wildcards):
	#print(expand("raw-data/{x}",x=config["samples"][wildcards.sample_name][wildcards.direction] ))
	return(expand("{x}",x=config["mrna_samples"][wildcards.mrna_sample]["r2"] ))


for i in 

rule all:
	input:
		["analysis/3-generate-fragment-lib/barcode-allele.tsv",
		"analysis/4-symmetric-barcode-quantification/startposcounts.lib.txt",
		"analysis/4-symmetric-barcode-quantification/startposcounts.481.txt"]

		
rule prep_merge_replicates:
	input:
		get_replicates
	output:
		temp("analysis/1-cluster-consensus/{sample_name}.{direction}.fastq")
	shell:
		"""
		zcat {input} > {output}
		"""
rule trim_for_calib:
	input:
		"analysis/1-cluster-consensus/{sample_name}.{direction}.fastq"
	output:
		temp("analysis/1-cluster-consensus/trimmed.{sample_name}.{direction}.fastq")
	shell:
		"""
		seqtk trimfq -L 75 {input} > {output}
		"""

rule calib_cluster:
	input:
		r1="analysis/1-cluster-consensus/trimmed.{sample_name}.r1.fastq",
		r2="analysis/1-cluster-consensus/trimmed.{sample_name}.r2.fastq"
	output:
		"analysis/1-cluster-consensus/{sample_name}.cluster"
	params:
		"-e 1 -k 5 -m 6 -l1 3 -l2 3 -t 2 -c 4 --no-sort --silent"
	shell:
		"""
		calib/calib -f {input.r1} -r {input.r2} -o analysis/{wildcards.sample_name}. {params}
		"""

rule calib_cons:
	input:
		r1="analysis/1-cluster-consensus/{sample_name}.r1.fastq",
		r2="analysis/1-cluster-consensus/{sample_name}.r2.fastq",
		cluster="analysis/1-cluster-consensus/{sample_name}.cluster"
	output:
		r1_cons_tmp=temp("analysis/1-cluster-consensus/constemp.{sample_name}.r1.fastq"),
		r2_cons_tmp=temp("analysis/1-cluster-consensus/constemp.{sample_name}.r2.fastq"),
		r1_cons_msa_tmp=temp("analysis/1-cluster-consensus/constemp.{sample_name}.r1.msa"),
		r2_cons_msa_tmp=temp("analysis/1-cluster-consensus/constemp.{sample_name}.r2.msa")
	shell:
		"""
		calib/consensus/calib_cons -q {input.r1} -q {input.r2} -o analysis/constemp.{wildcards.sample_name}.r1 analysis/constemp.{wildcards.sample_name}.r2 -c {input.cluster}
		"""

rule label_reads:
	input:
		r1="analysis/1-cluster-consensus/constemp.{sample_name}.{direction}.fastq",
		#r2="analysis/constemp.{sample_name}.r2.fastq"
	output:
		r1="analysis/1-cluster-consensus/cons.{sample_name}.{direction}.fastq",
		#r2="analysis/cons.{sample_name}.r2.fastq"
	shell:
		"""
		cat {input.r1}  | awk 'BEGIN {{OFS="\t"}} {{if(NR % 4 == 1) {{print $1_\"{wildcards.sample_name}\",$2;}} else print $0 }};' > {output.r1}
		"""


rule gather_reads:
	input:
		a="analysis/1-cluster-consensus/cons.longshort.{direction}.fastq",
		b="analysis/1-cluster-consensus/cons.shortlong.{direction}.fastq"
		#expand("analysis/cons.{sample_name}.{direction}.fastq", sample={wildcards.sample_name},direction={wildcards.direction})
		#r1=lambda wildcards: ["analysis/cons.{1}.r1.fastq".format(wildcards.sample_name)]
	output:
		r1="analysis/1-cluster-consensus/merged.{direction}.fastq"
	shell:
		"""
		cat {input.a} {input.b}  > {output.r1}
		"""

rule match_reads:
	input:
		r1="analysis/1-cluster-consensus/merged.r1.fastq",
		r2="analysis/1-cluster-consensus/merged.r2.fastq"
	output:
		orphan="analysis/2-match-reads/orphan.fastq",
		clustered="analysis/2-match-reads/clustered.interleaved.fastq",
		problematic_samesame="analysis/2-match-reads/problematic_samesame.interleaved.fastq",
		problematic_multiple="analysis/2-match-reads/problematic_multiple.interleaved.fastq",
		master_file="analysis/2-match-reads/master-barcode-cid.txt"
	shell:
		"""
		python code/snp-starrseq/read_id_matcher.py --f {input.r1} --r {input.r2} --p analysis/match-reads/
		"""

rule prep_collapse_fragments:
	input:
		"analysis/2-match-reads/clustered.interleaved.fastq"
	output:
		clust_r1=temp("analysis/3-generate-fragment-lib/clustered.r1.fastq"),
		clust_r2=temp("analysis/3-generate-fragment-lib/clustered.r2.fastq"),
	conda: 
		"env/env.collapse.yaml"
	shell:
		"""
		reformat.sh -Xmx100g in1={input} out1={output.clust_r1} out2={output.clust_r2} 
		"""

rule pre_collapse_trim:
	input:
		clust_r1="analysis/3-generate-fragment-lib/clustered.r1.fastq",
		clust_r2="analysis/3-generate-fragment-lib/clustered.r2.fastq"
	output:
		trim_r1="analysis/3-generate-fragment-lib/trimmed.clustered.r1.fastq",
		trim_r2="analysis/3-generate-fragment-lib/trimmed.clustered.r2.fastq",
	shell:
		"""
			seqtk trimfq -q 0.05 {input.clust_r1} > {output.trim_r1}
			seqtk trimfq -q 0.05 {input.clust_r2} > {output.trim_r2}
		"""


rule collapse_fragments:
	input:
		trim_r1="analysis/3-generate-fragment-lib/trimmed.clustered.r1.fastq",
		trim_r2="analysis/3-generate-fragment-lib/trimmed.clustered.r2.fastq",
	output:
		collapsed="analysis/3-generate-fragment-lib/collapsed-fragments.fastq"
	conda: 
		"env/env.collapse.yaml"
	shell:
		"""
		bbmerge.sh -Xmx100G ecco=t merge=t in1={input.trim_r1} in2={input.trim_r2} out={output.collapsed} 
		"""

rule map_collapsed_to_ref:
	input: 
		"analysis/3-generate-fragment-lib/collapsed-fragments.fastq"
	output:
		"analysis/3-generate-fragment-lib/collapsed-fragments.bam"
	conda:
		"env/env.mapping.yaml"
	shell:
		"""
		bwa mem -t 10 ~/tools/bwa_hg19_index/hg19.fa {input} | samtools view -Sb | samtools sort > {output} && samtools index {output}
		"""

rule create_variant_table:
	input:
		"analysis/3-generate-fragment-lib/collapsed-fragments.bam"
	output:
		"analysis/3-generate-fragment-lib/barcode-variant-table.tsv"
	conda:
		"env/env.variant-db.yaml"
	shell:
		"""
		python code/snp-starrseq/create-barcode-variant-db.py --bam_input {input} > {output}
		"""

rule create_mutation_database:
	input:
		bam="analysis/3-generate-fragment-lib/collapsed-fragments.bam",
		variant_table="analysis/3-generate-fragment-lib/barcode-variant-table.tsv"
	output:
		"analysis/3-generate-fragment-lib/barcode-allele.tsv"
	conda:
		"env/env.variant-db.yaml"
	params:
		bed=config["capture_bed"]
	shell:
		"""
		python code/snp-starrseq/dev-segment-variant-matrix.py --bam_input {input.bam} --var_db {input.variant_table} --csv_out {output} --mat_out False --bed_file {params.bed}
		"""
rule mrna_quantification_startpos:
	input:
		r1=get_mrna_r1,
		r2=get_mrna_r2,
	output:
		countx="analysis/4-symmetric-barcode-quantification/startposcounts.{mrna_sample}.txt",
		bam="analysis/4-symmetric-barcode-quantification/{mrna_sample}.bam",
		sortedbam="analysis/4-symmetric-barcode-quantification/{mrna_sample}.sorted.bam"
	shell:
		"""
		bwa mem -t60 ~/tools/bwa_hg19_index/hg19.fa {input.r1} {input.r2} | samtools view -Sb > {output.bam}
		samtools view -f2 -F2304 {output.bam} | python code/snp-starrseq/create-umi-directory-paired.py - | cut -f 1,2 | sort -T . --parallel=10 | uniq -c | awk '{{print $1,$2,$3}}' > {output.countx}
		samtools sort -@10 -m10G {output.bam} > {output.sortedbam}
		samtools index {output.sortedbam}
		"""


rule longread_quantification_startpos:
	input:
		config["longread"]
	output:
		"analysis/5-longread/longread-counts.txt"
	shell:
		"""
		samtools view -q10 -F2052 {input} | python code/snp-starrseq/create-umi-directory-longread.py | cut -f 1,2 | sort | uniq -c  | awk '{print $1,$2,$3}' > {output}
		"""