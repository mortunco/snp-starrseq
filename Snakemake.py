def get_replicates(wildcards):
	#print(expand("raw-data/{x}",x=config["samples"][wildcards.sample_name][wildcards.direction] ))
	return(expand("{x}",x=config["dna_samples"][wildcards.sample_name][wildcards.direction] ))
def get_mrna_r1(wildcards):
	#print(expand("raw-data/{x}",x=config["samples"][wildcards.sample_name][wildcards.direction] ))
	return(expand("{x}",x=config["mrna_samples"][wildcards.mrna_sample]["r1"] ))
def get_mrna_r2(wildcards):
	#print(expand("raw-data/{x}",x=config["samples"][wildcards.sample_name][wildcards.direction] ))
	return(expand("{x}",x=config["mrna_samples"][wildcards.mrna_sample]["r2"] ))


# rule all:
# 	input:
# 		["analysis/collapsed/barcode-allele.tsv"]

rule all:
	input:
		["analysis/collapsed/barcode-allele.tsv",
		"analysis/mrna/barcodecounts.lib.txt",
		"analysis/mrna/barcodecounts.481.txt",
		"analysis/mrna/startposcounts.lib.txt",
		"analysis/mrna/startposcounts.481.txt"]

		
rule prep_merge_replicates:
	input:
		get_replicates
	output:
		temp("analysis/{sample_name}.{direction}.fastq")
	shell:
		"""
		zcat {input} > {output}
		"""
rule trim_for_calib:
	input:
		"analysis/{sample_name}.{direction}.fastq"
	output:
		temp("analysis/trimmed.{sample_name}.{direction}.fastq")
	shell:
		"""
		seqtk trimfq -L 75 {input} > {output}
		"""

rule calib_cluster:
	input:
		r1="analysis/trimmed.{sample_name}.r1.fastq",
		r2="analysis/trimmed.{sample_name}.r2.fastq"
	output:
		"analysis/{sample_name}.cluster"
	params:
		"-e 1 -k 5 -m 6 -l1 3 -l2 3 -t 2 -c 4 --no-sort --silent"
	shell:
		"""
		calib/calib -f {input.r1} -r {input.r2} -o analysis/{wildcards.sample_name}. {params}
		"""

rule calib_cons:
	input:
		r1="analysis/{sample_name}.r1.fastq",
		r2="analysis/{sample_name}.r2.fastq",
		cluster="analysis/{sample_name}.cluster"
	output:
		r1_cons_tmp=temp("analysis/constemp.{sample_name}.r1.fastq"),
		r2_cons_tmp=temp("analysis/constemp.{sample_name}.r2.fastq"),
		r1_cons_msa_tmp=temp("analysis/constemp.{sample_name}.r1.msa"),
		r2_cons_msa_tmp=temp("analysis/constemp.{sample_name}.r2.msa")
	shell:
		"""
		calib/consensus/calib_cons -q {input.r1} -q {input.r2} -o analysis/constemp.{wildcards.sample_name}.r1 analysis/constemp.{wildcards.sample_name}.r2 -c {input.cluster}
		"""

rule label_reads:
	input:
		r1="analysis/constemp.{sample_name}.{direction}.fastq",
		#r2="analysis/constemp.{sample_name}.r2.fastq"
	output:
		r1="analysis/cons.{sample_name}.{direction}.fastq",
		#r2="analysis/cons.{sample_name}.r2.fastq"
	shell:
		"""
		cat {input.r1}  | awk 'BEGIN {{OFS="\t"}} {{if(NR % 4 == 1) {{print $1_\"{wildcards.sample_name}\",$2;}} else print $0 }};' > {output.r1}
		"""


rule gather_trimquality_reads:
	input:
		a="analysis/cons.longshort.{direction}.fastq",
		b="analysis/cons.shortlong.{direction}.fastq"
		#expand("analysis/cons.{sample_name}.{direction}.fastq", sample={wildcards.sample_name},direction={wildcards.direction})
		#r1=lambda wildcards: ["analysis/cons.{1}.r1.fastq".format(wildcards.sample_name)]
	output:
		r1="analysis/merged.{direction}.fastq"
	shell:
		"""
		cat {input.a} {input.b}  > {output.r1}
		"""

rule match_reads:
	input:
		r1="analysis/merged.r1.fastq",
		r2="analysis/merged.r2.fastq"
	output:
		orphan="analysis/match-reads/orphan.fastq",
		clustered="analysis/match-reads/clustered.interleaved.fastq",
		problematic_samesame="analysis/match-reads/problematic_samesame.interleaved.fastq",
		problematic_multiple="analysis/match-reads/problematic_multiple.interleaved.fastq",
		master_file="analysis/match-reads/master-barcode-cid.txt"
	shell:
		"""
		python code/snp-starrseq/read_id_matcher.py --f {input.r1} --r {input.r2} --p analysis/match-reads/
		"""

rule collapse_fragments:
	input:
		"analysis/match-reads/clustered.interleaved.fastq"
	output:
		clust_r1=temp("analysis/collapsed/clustered.r1.fastq"),
		clust_r2=temp("analysis/collapsed/clustered.r2.fastq"),
	conda: 
		"env/env.collapse.yaml"
	shell:
		"""
		reformat.sh -Xmx100g in1={input} out1={output.clust_r1} out2={output.clust_r2} 
		"""

rule post_collapse_trim:
	input:
		clust_r1="analysis/collapsed/clustered.r1.fastq",
		clust_r2="analysis/collapsed/clustered.r2.fastq"
	output:
		trim_r1="analysis/collapsed/trimmed.clustered.r1.fastq",
		trim_r2="analysis/collapsed/trimmed.clustered.r2.fastq",
	shell:
		"""
			seqtk trimfq -q 0.05 {input.clust_r1} > {output.trim_r1}
			seqtk trimfq -q 0.05 {input.clust_r2} > {output.trim_r2}
		"""


rule collapse_fragments_prep:
	input:
		trim_r1="analysis/collapsed/trimmed.clustered.r1.fastq",
		trim_r2="analysis/collapsed/trimmed.clustered.r2.fastq",
	output:
		collapsed="analysis/collapsed/collapsed-fragments.fastq"
	conda: 
		"env/env.collapse.yaml"
	shell:
		"""
		bbmerge.sh -Xmx100G ecco=t merge=t in1={input.trim_r1} in2={input.trim_r2} out={output.collapsed} 
		"""

rule mapping_collapsed_ref:
	input: 
		"analysis/collapsed/collapsed-fragments.fastq"
	output:
		"analysis/collapsed/collapsed.bam"
	conda:
		"env/env.mapping.yaml"
	shell:
		"""
		bwa mem -t 10 ~/tools/bwa_hg19_index/hg19.fa {input} | samtools view -Sb | samtools sort > {output} && samtools index {output}
		"""

rule create_variant_table:
	input:
		"analysis/collapsed/collapsed.bam"
	output:
		"analysis/collapsed/barcode-variant-table.tsv"
	conda:
		"env/env.variant-db.yaml"
	shell:
		"""
		python code/snp-starrseq/create-barcode-variant-db.py --bam_input {input} > {output}
		"""

rule create_mutation_database:
	input:
		bam="analysis/collapsed/collapsed.bam",
		variant_table="analysis/collapsed/barcode-variant-table.tsv"
	output:
		"analysis/collapsed/barcode-allele.tsv"
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
		countx="analysis/mrna/startposcounts.{mrna_sample}.txt",
		bam="analysis/mrna/{mrna_sample}.bam",
		sortedbam="analysis/mrna/{mrna_sample}.sorted.bam"
	shell:
		"""
		bwa mem -t60 ~/tools/bwa_hg19_index/hg19.fa {input.r1} {input.r2} | samtools view -Sb > {output.bam}
		samtools view -f2 -F2304 {output.bam} | python code/snp-starrseq/create-umi-directory-paired.py - | cut -f 1,2 | sort -T . --parallel=10 | uniq -c | awk '{{print $1,$2,$3}}' > {output.countx}
		samtools sort -@10 -m10G {output.bam} > {output.sortedbam}
		samtools index {output.sortedbam}
		"""

		

rule mrna_quantification_barcode:
	input:
		r1=get_mrna_r1,
		r2=get_mrna_r2,
	output:
		"analysis/mrna/barcodecounts.{mrna_sample}.txt"
	conda:
		"env/env.mapping.yaml"
	shell:
		"""
		paste -d "\t"  \
		<(zcat {input.r1} | paste -d "\t" - - - - ) \
		<(zcat {input.r2} | paste -d "\t" - - - - )  | cut -f  1,2,6 | awk '{{print substr($2,0,12) substr($3,0,12)}}' | sort | uniq -c | awk '{{print $2,$1}}' | sort -k2,2n  > {output}
		"""
