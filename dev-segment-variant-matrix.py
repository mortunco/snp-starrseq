#!/usr/bin/env python
# coding: utf-8




import pysam
import numpy
import pandas
import argparse



# required arg

parser = argparse.ArgumentParser()
   
parser.add_argument('--bam_input', required=True)
parser.add_argument("--var_db", required=True)
parser.add_argument('--csv_out', required=True)
parser.add_argument('--mat_out', required=True)
parser.add_argument('--bed_file', required=True)

args = parser.parse_args()





db_file=args.var_db
bam_file=args.bam_input
output=args.csv_out
matrix_output=args.mat_out
target_region=args.bed_file




# db_file="/groups/lackgrp/projects/inv-mattfreedman-snpstarrseq/analysis/process-asymmetric/mutation-database-v2/segment-chr1_10268688_10273174-db.txt"
# bam_file="/groups/lackgrp/projects/inv-mattfreedman-snpstarrseq/analysis/process-asymmetric/mutation-database-v2/segment-chr1_10268688_10273174.bam"
# output="/groups/lackgrp/projects/inv-mattfreedman-snpstarrseq/analysis/process-asymmetric/mutation-database-v2/segment-chr10_104417350_104419349-matrix.tsv"




print("Creating Dictionaries done")
mutation_dict={} ### keeps all alleles of same position mutations 
mutation_position=[] ### keeps positions (mutations) to be used by pile up ###
barcode_index_dict={} ### keeps collapsed fragent barcode --> index information.
with open(db_file, 'r') as file_in:
    ### temp example ['ATGTGTAAAATATAAACAAAACAT', 'chr6:109326943-ATGCAT', 'ATGTTTTGTTTATATTTTACACAT', 'snp', 'chr6', '109326945', 'C', 'G', 'PASS']
    for line in file_in:
        temp=line.rstrip("\n").split("\t")
        index, fragment_name, event_type, event_chr, event_pos, event_ref, event_alt, temp_x = temp

        if fragment_name not in barcode_index_dict:
            barcode_index_dict[fragment_name]=index
        if event_type != "snp" or event_alt == "N": ### I discarded all the indel events now  + REMOVES all N cases ###
            continue
        if event_chr + ";" + event_pos not in mutation_dict:
            mutation_dict[event_chr + ";" + event_pos] = [event_ref+";"+event_alt]
        else:
            if mutation_dict[event_chr + ";" + event_pos] != [event_ref+";"+event_alt]:
                mutation_dict[event_chr + ";" + event_pos].append(event_ref+";"+event_alt)
            else:
                pass
        #mutation_position.append(event_chr + ";" + event_pos)
        #mutation_list.append(";".join(temp[4:8]))
        #mutation_list.append(";".join([event_chr,event_pos,event_ref,event_alt]))

#{'chrX;66765158': ['TGCAGCAGCA;T', 'T;TGCA', 'TGCAGCA;T', 'TGCAGCAGCAGCA;T', 'T;TGCAGCA', 'T;TGCAGCAGCA', 'TGCA;T', 'TGCAGCAGCAGCAGCAGCAGCA;T', 'TGCAGCAGCAGCAGCAGCAGCAGCA;T', 'T;TGCAGCAGCAGCA']
blacklist_mutations=[]
for event_position in list(mutation_dict.keys()):
    if len(mutation_dict[event_position]) >=3:
        mutation_dict.pop(event_position,None)

# for i in blacklist_mutations:
#     mutation_dict.pop(i,None)
#mutation_position = [i for i in mutation_position if i not in blacklist_mutations]
mutation_position= list(mutation_dict.keys())
print(mutation_dict["chr6;109326183"])
print("Filtering Blacklist done")
mydf=list()
with open(target_region, 'r') as file_in:
    segment_no=1
    ### for example
    ### chr1    10438799        10440798        1       snp_rs41310365
    for bed_line in file_in:
        capture_chr,capture_start,capture_end,temp_x,temp_y = bed_line.rstrip("\n").split("\t")
        segment=[capture_chr,capture_start,capture_end] # chr start end of bed
        print("Analysing Segment no {}      ".format(segment_no,"_".join(segment)))
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for pileupcolumn in samfile.pileup(segment[0], int(segment[1]),int(segment[2])):
            if pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1) in mutation_position:#["chr10;104418948","chr10;104418946","chr10;104418945"]:     
                for allele in mutation_dict[pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1)]: ### we are doing this because some events are multi alallic therefore i have to go through all of them.
                    #print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                    for pileupread in pileupcolumn.pileups:
                        # dont delete this this is required for inspection if stuff goes bad
                        # print("allele : {} | qname: {} | qpos : {} | pos : {} | indel : {} | is_del : {} | is_refskip : {} ".format(allele,pileupread.alignment.query_name,
                        #        pileupread.query_position,
                        #        pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1),
                        #        pileupread.indel,
                        #        pileupread.is_del,
                        #        pileupread.is_refskip)
                        #       )                        
                        if pileupread.indel != 0: 
                        ## Check if the mutation we are genotyping if indel, if this base is indel, then it directly supports indels.
                        ## I used this line in the original version that I also calculated the genotype of the indel. But that is complex. Therefore, its depracated.
                        ## I am commanding out below, because this would have added any read with indel to the support of SNPs. which is not correct
                            pass
                            #temp=[barcode_index_dict[pileupread.alignment.qname],pileupread.alignment.qname,pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1) + ";" + allele, 2]
                            #mydf.append(temp)
                        else:
                            if pileupread.query_position == None: 
                            ## This is the case, when indel starts but because its deletion the sequence is none in the read. Thefore, I am not discarding this case.
                                pass #### 
                            else:
                                ### If no indel, then it checks the first base, if mutation was indel checks if it matches witht he first reference. it should be same.
                                if pileupread.alignment.query_sequence[pileupread.query_position] in allele.split(";")[0]:
                                    temp=[barcode_index_dict[pileupread.alignment.qname],pileupread.alignment.qname,pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1) + ";" + allele,1]
                                    mydf.append(temp)
                                else:
                                    temp=[barcode_index_dict[pileupread.alignment.qname],pileupread.alignment.qname,pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1) + ";" + allele,2]
                                    #print("\t"+" ".join(temp))
                                    mydf.append(temp)
        segment_no+=1                
        samfile.close()

x=pandas.DataFrame(mydf,columns =["start_aln_UMI",'fragment_name', "mutation","type"])
if matrix_output == True:
    x.pivot_table(index=['fragment_name'], columns='mutation',values='type').to_csv(output,sep="\t")
else:
    x.to_csv(output,sep="\t", index=False)
print("Done")

