#!/usr/bin/env python
# coding: utf-8




import pysam
import numpy
import pandas
import argparse





import argparse

# required arg

parser = argparse.ArgumentParser()
   
parser.add_argument('--bam_input', required=True)
parser.add_argument("--var_db", required=True)
parser.add_argument('--csv_out', required=True)
parser.add_argument('--mat_out', required=True)

args = parser.parse_args()





db_file=args.var_db
bam_file=args.bam_input
output=args.csv_out
matrix_output=args.mat_out





# db_file="/groups/lackgrp/projects/inv-mattfreedman-snpstarrseq/analysis/process-asymmetric/mutation-database-v2/segment-chr1_10268688_10273174-db.txt"
# bam_file="/groups/lackgrp/projects/inv-mattfreedman-snpstarrseq/analysis/process-asymmetric/mutation-database-v2/segment-chr1_10268688_10273174.bam"
# output="/groups/lackgrp/projects/inv-mattfreedman-snpstarrseq/analysis/process-asymmetric/mutation-database-v2/segment-chr10_104417350_104419349-matrix.tsv"





mutation_dict={} ### keeps all alleles of same position mutations 
mutation_position=[] ### keeps positions (mutations) to be used by pile up ###
mutation_list=[] ### keeps tracks of mutations so that we dont have duplicates of same mutations from different barcodes 
with open(db_file, 'r') as file_in:
    for line in file_in:
        temp=line.rstrip("\n").split("\t")
        if temp[1] == "WT":
            continue
        #print(temp)
        if ";".join(temp[2:4]) not in mutation_dict:
                mutation_dict[";".join(temp[2:4])] = [temp[4]+";"+temp[5]]
        elif ";".join(temp[2:6]) in mutation_list:
            continue
        else:
            mutation_dict[";".join(temp[2:4])].append(temp[4]+";"+temp[5])
        mutation_position.append(";".join(temp[2:4]))
        mutation_list.append(";".join(temp[2:6]))





blacklist_mutations=[]
for i,j in mutation_dict.items():
    if len(j) >=3:
        blacklist_mutations.append(i)





for i in blacklist_mutations:
    mutation_dict.pop(i,None)
mutation_position = [i for i in mutation_position if i not in blacklist_mutations]
mutation_list = [i for i in mutation_position if i not in mutation_list]





mydf=list()
samfile = pysam.AlignmentFile(bam_file, "rb")
for pileupcolumn in samfile.pileup():
    if pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1) in mutation_position:#["chr10;104418948","chr10;104418946","chr10;104418945"]:     
        for allele in mutation_dict[pileupcolumn.reference_name + ";" + str(pileupcolumn.pos+1)]: ### we are doing this becausem some events are multi alallic therefore i have to go through all of them.
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



                
samfile.close()

x=pandas.DataFrame(mydf,columns =['barcode', "mutation","type"])
if matrix_output == True:
    x.pivot_table(index=['barcode'], columns='mutation',values='type').to_csv(output,sep="\t")
else:
    x.to_csv(output,sep="\t", index=False)
print("done xd")

