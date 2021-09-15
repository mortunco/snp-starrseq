def get_replicates(wildcards):
	return(expand("{x}",x=config["asym_samples"][wildcards.sample_name][wildcards.direction] ))
def get_mrna_r1(wildcards):
	return(expand("{x}",x=config["symmetric_samples"][wildcards.mrna_sample]["r1"] ))
def get_mrna_r2(wildcards):
	return(expand("{x}",x=config["symmetric_samples"][wildcards.mrna_sample]["r2"] ))

def genereate_result():
	results=[]
	run_name=config["run_name"]
	if config["symmetric_samples"] != "":
		for i in list(config["symmetric_samples"].keys()):
			results.append("analysis/{}/4-symmetric-barcode-quantification/startposcounts.{}.txt".format(run_name,i))
	if config["longread_samples"] != "" and "longread_samples" in config:
		results.append("analysis/{}/5-longread/longread-counts.txt".format(run_name))
	results.append("analysis/{}/3-generate-fragment-lib/barcode-allele.tsv".format(run_name))
	return results

rule all:
	input: genereate_result()

rule prep_merge_replicates:
	input:
		get_replicates
	output:
		temp("analysis/{run_name}/1-cluster-consensus/{sample_name}.{direction}.fastq")
	conda:
		"env/env.seqtk.yaml"
	threads:1
	shell:
		"""
		zcat {input} > {output}
		"""
rule trim_for_calib:
	input:
		"analysis/{run_name}/1-cluster-consensus/{sample_name}.{direction}.fastq"
	output:
		temp("analysis/{run_name}/1-cluster-consensus/trimmed.{sample_name}.{direction}.fastq")
	conda:
		"env/env.seqtk.yaml"
	threads:1
	shell:
		"""
		seqtk trimfq -L 75 {input} > {output}
		"""

rule calib_cluster:
	input:
		r1="analysis/{run_name}/1-cluster-consensus/trimmed.{sample_name}.r1.fastq",
		r2="analysis/{run_name}/1-cluster-consensus/trimmed.{sample_name}.r2.fastq"
	output:
		"analysis/{run_name}/1-cluster-consensus/{sample_name}.cluster"
	params:
		config["calib_params"]
	log: 
		"logs/{run_name}/1-calib-cluster.{sample_name}.log"
	threads:4
	shell:
		"""
		./code/calib/calib/calib -f {input.r1} -r {input.r2} -o analysis/{wildcards.run_name}/1-cluster-consensus/{wildcards.sample_name}. {params} 2> {log}
		"""

rule calib_cons:
	input:
		r1="analysis/{run_name}/1-cluster-consensus/{sample_name}.r1.fastq",
		r2="analysis/{run_name}/1-cluster-consensus/{sample_name}.r2.fastq",
		cluster="analysis/{run_name}/1-cluster-consensus/{sample_name}.cluster"
	output:
		r1_cons_tmp=temp("analysis/{run_name}/1-cluster-consensus/constemp.{sample_name}.r1.fastq"),
		r2_cons_tmp=temp("analysis/{run_name}/1-cluster-consensus/constemp.{sample_name}.r2.fastq"),
		r1_cons_msa_tmp=temp("analysis/{run_name}/1-cluster-consensus/constemp.{sample_name}.r1.msa"),
		r2_cons_msa_tmp=temp("analysis/{run_name}/1-cluster-consensus/constemp.{sample_name}.r2.msa")
	log: 
		"logs/{run_name}/2-calib-cons.{sample_name}.log"
	threads: 4
	shell:
		"""
		./code/calib/calib/consensus/calib_cons -q {input.r1} -q {input.r2} -o analysis/{wildcards.run_name}/1-cluster-consensus/constemp.{wildcards.sample_name}.r1 analysis/{wildcards.run_name}/1-cluster-consensus/constemp.{wildcards.sample_name}.r2 -c {input.cluster} 2> {log}
		"""

rule label_reads:
	input:
		r1="analysis/{run_name}/1-cluster-consensus/constemp.{sample_name}.{direction}.fastq",
	output:
		r1="analysis/{run_name}/1-cluster-consensus/cons.{sample_name}.{direction}.fastq",
	threads:1
	shell:
		"""
		cat {input.r1}  | awk 'BEGIN {{OFS="\t"}} {{if(NR % 4 == 1) {{print $1_\"{wildcards.sample_name}\",$2;}} else print $0 }};' > {output.r1}
		"""


rule gather_reads:
	input:
		a="analysis/{run_name}/1-cluster-consensus/cons.longshort.{direction}.fastq",
		b="analysis/{run_name}/1-cluster-consensus/cons.shortlong.{direction}.fastq"
	output:
		r1="analysis/{run_name}/1-cluster-consensus/merged.{direction}.fastq"
	threads:1
	shell:
		"""
		cat {input.a} {input.b}  > {output.r1}
		"""

rule match_reads:
	input:
		r1="analysis/{run_name}/1-cluster-consensus/merged.r1.fastq",
		r2="analysis/{run_name}/1-cluster-consensus/merged.r2.fastq"
	output:
		orphan="analysis/{run_name}/2-match-reads/orphan.fastq",
		clustered_r1="analysis/{run_name}/2-match-reads/clustered.r1.fastq",
		clustered_r2="analysis/{run_name}/2-match-reads/clustered.r2.fastq",
		problematic_samesame="analysis/{run_name}/2-match-reads/problematic_samesame.interleaved.fastq",
		problematic_multiple="analysis/{run_name}/2-match-reads/problematic_multiple.interleaved.fastq",
		master_file="analysis/{run_name}/2-match-reads/master-barcode-cid.txt"
	threads:1
	shell:
		"""
		python code/snp-starrseq/read_id_matcher.py --f {input.r1} --r {input.r2} --p analysis/{wildcards.run_name}/2-match-reads/
		"""

rule pre_collapse_trim:
	input:
		clust_r1="analysis/{run_name}/2-match-reads/clustered.r1.fastq",
		clust_r2="analysis/{run_name}/2-match-reads/clustered.r2.fastq"
	output:
		trim_r1="analysis/{run_name}/3-generate-fragment-lib/trimmed.clustered.r1.fastq",
		trim_r2="analysis/{run_name}/3-generate-fragment-lib/trimmed.clustered.r2.fastq",
	conda:
		"env/env.seqtk.yaml"
	threads:1
	shell:
		"""
			seqtk trimfq -q 0.05 {input.clust_r1} > {output.trim_r1}
			seqtk trimfq -q 0.05 {input.clust_r2} > {output.trim_r2}
		"""


rule collapse_fragments:
	input:
		trim_r1="analysis/{run_name}/3-generate-fragment-lib/trimmed.clustered.r1.fastq",
		trim_r2="analysis/{run_name}/3-generate-fragment-lib/trimmed.clustered.r2.fastq",
	output:
		collapsed=temp("analysis/{run_name}/3-generate-fragment-lib/collapsed-fragments.fastq")
	conda: 
		"env/env.collapse.yaml"
	log: 
		"logs/{run_name}/3-bbmerge.log"
	threads:1
	shell:
		"""
		bbmerge.sh -Xmx100G ecco=t merge=t in1={input.trim_r1} in2={input.trim_r2} out={output.collapsed} 2>{log}
		"""

rule map_collapsed_to_ref:
	input: 
		"analysis/{run_name}/3-generate-fragment-lib/collapsed-fragments.fastq"
	output:
		"analysis/{run_name}/3-generate-fragment-lib/collapsed-fragments.bam"
	conda:
		"env/env.mapping.yaml"
	log: 
		"logs/{run_name}/3-map-fragments.log"
	params:
		config["bwa_ref"]
	threads:16
	shell:
		"""
		bwa mem -t{threads} {params} {input} | samtools view -Sb | samtools sort -@{threads} -m10G> {output} && samtools index {output} 2> {log}
		"""

rule create_variant_table:
	input:
		"analysis/{run_name}/3-generate-fragment-lib/collapsed-fragments.bam"
	output:
		"analysis/{run_name}/3-generate-fragment-lib/barcode-variant-table.tsv"
	conda:
		"env/env.variant-db.yaml"
	threads:1
	shell:
		"""
		python code/snp-starrseq/create-barcode-variant-db.py --bam_input {input} > {output}
		"""

rule create_mutation_database:
	input:
		bam="analysis/{run_name}/3-generate-fragment-lib/collapsed-fragments.bam",
		variant_table="analysis/{run_name}/3-generate-fragment-lib/barcode-variant-table.tsv"
	output:
		"analysis/{run_name}/3-generate-fragment-lib/barcode-allele.tsv"
	conda:
		"env/env.variant-db.yaml"
	params:
		bed=config["capture_bed"]
	log: 
		"logs/{run_name}/3-mutation-database.log"
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
		sortedbam="analysis/{run_name}/4-symmetric-barcode-quantification/{mrna_sample}.sorted.bam"
	conda:
		"env/env.mapping.yaml"	
	params:
		config["bwa_ref"]
	threads:16
	shell:
		"""
		bwa mem -t{threads} {params} {input.r1} {input.r2} | samtools view -Sb > {output.bam}
		samtools view -f2 -F2048 {output.bam} | python code/snp-starrseq/create-umi-directory-paired.py - | cut -f 1 | sort -T . --parallel={threads} | uniq -c | awk '{{print $1,$2}}' > {output.countx}
		samtools sort -@{threads} -m10G {output.bam} > {output.sortedbam}
		samtools index {output.sortedbam}
		"""
# rule NB
# shell:
# 	"""
# 	python3 allele_specific_enhancer_activity.py \
#     --dna=yi-results-new/startposcounts.SNP_CL.txt \
#     --rna=yi-results-new/startposcounts.SNPR1_48.txt,yi-results-new/startposcounts.SNPR2_48.txt,yi-results-new/startposcounts.SNPR3_48.txt \
#     --barcode_allele=yi-results-new/barcode-allele.tsv --output_prefix=48h 
# 	"""

rule longread_quantification_startpos:
	input:
		config["longread_samples"]
	output:
		"analysis/{run_name}/5-longread/longread-counts.txt"
	conda:
		"env/env.mapping.yaml"	
	threads: 16
	shell:
		"""
		samtools view -q2 -F2052 {input} | python code/snp-starrseq/create-umi-directory-longread.py | cut -f 1 | sort -T . --parallel={threads} | uniq -c  | awk '{{print $1,$2}}' > {output}
		"""
