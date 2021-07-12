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




print("Creating Dictionaries")
mutation_dict={} ### keeps all alleles of same position mutations 
mutation_position=[] ### keeps positions (mutations) to be used by pile up ###
mutation_list=[] ### keeps tracks of mutations so that we dont have duplicates of same mutations from different barcodes 
with open(db_file, 'r') as file_in:
    ### temp example ['ATGTGTAAAATATAAACAAAACAT', 'chr6:109326943-ATGCAT', 'ATGTTTTGTTTATATTTTACACAT', 'snp', 'chr6', '109326945', 'C', 'G', 'PASS']
    for line in file_in:
        temp=line.rstrip("\n").split("\t")
        #print(temp)
        if temp[3] != "snp":
            continue
        if ";".join(temp[3:5]) not in mutation_dict:
            mutation_dict[";".join(temp[4:6])] = [temp[6]+";"+temp[7]]
        elif ";".join(temp[3:7]) in mutation_list:
            continue
        else:
            mutation_dict[";".join(temp[4:6])] = [temp[6]+";"+temp[7]]
        mutation_position.append(";".join(temp[4:6]))
        mutation_list.append(";".join(temp[4:8]))
print("Done")


#{'chrX;66765158': ['TGCAGCAGCA;T', 'T;TGCA', 'TGCAGCA;T', 'TGCAGCAGCAGCA;T', 'T;TGCAGCA', 'T;TGCAGCAGCA', 'TGCA;T', 'TGCAGCAGCAGCAGCAGCAGCA;T', 'TGCAGCAGCAGCAGCAGCAGCAGCA;T', 'T;TGCAGCAGCAGCA']
blacklist_mutations=[]
print("Filtering Blacklist")
for i,j in mutation_dict.items():
    if len(j) >=3:
        blacklist_mutations.append(i)

for i in blacklist_mutations:
    mutation_dict.pop(i,None)
mutation_position = [i for i in mutation_position if i not in blacklist_mutations]
mutation_list = [i for i in mutation_position if i not in mutation_list]
print("Done")
mydf=list()
with open(target_region, 'r') as file_in:
    segment_no=1
    ### for example
    ### chr1    10438799        10440798        1       snp_rs41310365
    for line in file_in:
        temp = line.rstrip("\n").split("\t")
        segment=[temp[0],temp[1],temp[2]] # chr start end of bed
        print("Analysing Segment no {}".format(segment_no,"_".join(segment)))
        samfile = pysam.AlignmentFile(bam_file, "rb")
        #print([segment[0], int(segment[1]),int(segment[2])])
        for pileupcolumn in samfile.pileup(segment[0], int(segment[1]),int(segment[2])):
            if pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1) in mutation_position:#["chr10;104418948","chr10;104418946","chr10;104418945"]:     
                for allele in mutation_dict[pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1)]: ### we are doing this because some events are multi alallic therefore i have to go through all of them.
                    #print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                    for pileupread in pileupcolumn.pileups:
        #                 print("allele : {} | qname: {} | pos : {} | seq : {}".format(allele,
        #                     pileupread.alignment.query_name,
        #                       pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1),
        #                       pileupread.alignment.query_sequence[pileupread.query_position])
        #                      )
        #                 print("qname: {} | qpos : {} | pos : {} | indel : {} | is_del : {} | is_refskip : {} ".format(pileupread.alignment.query_name,
        #                       pileupread.query_position,
        #                       pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1),
        #                       pileupread.indel,
        #                       pileupread.is_del,
        #                       pileupread.is_refskip)
        #                      )
                        
                        if pileupread.indel != 0: ## checks indel if its indel then its directly mutation.
                            temp=[pileupread.alignment.qname,pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1) + ";" + allele, 2]
                            #print("\t"+" ".join(temp))
                        else:
                            if pileupread.query_position == None:
                                pass #### 
                            else:
                                ### If no indel, then it checks the first base, if mutation was indel checks if it matches witht he first reference. it should be same.
                                if pileupread.alignment.query_sequence[pileupread.query_position] in allele.split(";")[0]:#mutation_ref[pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1)]:
                                    temp=[pileupread.alignment.qname,
                                          pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1) + ";" + allele,1]
                                    #print("\t"+" ".join(temp))
                                    mydf.append(temp)
                                    
                                else:
                                    temp=[pileupread.alignment.qname,pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1) + ";" + allele,2]
                                    #print("\t"+" ".join(temp))
                                    mydf.append(temp)



        segment_no+=1                
        samfile.close()

x=pandas.DataFrame(mydf,columns =['fragment_name', "mutation","type"])
if matrix_output == True:
    x.pivot_table(index=['fragment_name'], columns='mutation',values='type').to_csv(output,sep="\t")
else:
    x.to_csv(output,sep="\t", index=False)
print("done xd")

