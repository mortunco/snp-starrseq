import argparse
import os
import sys
import pandas as pd
import pybedtools as pybed


parser = argparse.ArgumentParser(description="Gets barcode of each snp event then separates two bam files (REF/ALT)")

parser.add_argument("-i","--input_vis_file",type=str,required=True,help="vis info file.")
parser.add_argument("-g","--genome",type=str,required=True,help="genome version")
parser.add_argument("-r","--result_table",type=str,required=True,help="bi-allelic result table")
parser.add_argument("-b","--barcode_allele",type=str,required=True,help="barcode allele file generated by either pacbio or asymmetricfile")
parser.add_argument("-o","--output_dir",type=str,help="output directory that contains all bedgraph files.")

def generate_bedgraph(index_event, event_type, snpid,chromsizes,event,output_dir,debug=False):
    if event_type == 1: ### REF allele
        reads=index_event[(index_event["mutation"] == event) & (index_event["type"] == event_type )][["chrid","start","end"]]
        reads=pybed.BedTool.from_dataframe(reads)
        reads=reads.genomecov(bg=True,g=chromsizes,scale=1/int(result_table[result_table["event"] == event]["ref_count"]))
        reads = reads.to_dataframe()
        if debug != True:
            reads.to_csv(f'{output_dir}/{snpid}_{event.replace(";","-")}.REF.bedGraph',sep="\t",header=False,index=False)
    elif event_type == 2: ### for VAR allele
        reads=index_event[(index_event["mutation"] == event) & (index_event["type"] == event_type )][["chrid","start","end"]]
        reads=pybed.BedTool.from_dataframe(reads)
        reads=reads.genomecov(bg=True,g=chromsizes,scale=1/int(result_table[result_table["event"] == event]["alt_count"]))
        reads = reads.to_dataframe() 
        if debug != True:
            reads.to_csv(f'{output_dir}/{snpid}_{event.replace(";","-")}.VAR.bedGraph',sep="\t",header=False,index=False)
    


args = parser.parse_args()
output_dir = args.output_dir
pybed.helpers.set_tempdir(output_dir)

# Barcode allale file containts fragment index to SNP. Type column is either 1 for WT , 2 for VAR 
# start_aln_UMI   fragment_name   mutation        type
# chr6:109321504-GGAATC   GGAGTAGTAATAACTGTAAAGATC        chr6;109321800;T;C      1
# chr6:109321541-GTAATC   GATCTCATAAGTTTACAGCACTAC        chr6;109321800;T;C      1
# chr6:109321555-CAACAC   GTGTCATTGACAACATTTTAGTTG        chr6;109321800;T;C      1
# chr6:109321556-GAACGT   ACGTTTAATCTGCACATTTAGTTC        chr6;109321800;T;C      1
index_event=args.barcode_allele
index_event=pd.read_table(index_event,sep="\t",header=0)
index_event=index_event.drop("fragment_name",axis=1)

### result table generated in the prevoius step.
# SNP     event   ref_count       alt_count       ref_expr        alt_expr        alt_effect      alt_freq        pvalue  z       padj
# rs10457185      chr6;109322477;G;C      46      8       2178    159     -0.7648405092458077     0.31846838022530427     0.23511709646737156 -1.1872802480887972     0.3134894619564954
# rs2273668       chr6;109323519;G;T      60      10      1998    111     -1.1019636195016487     0.2484278291072463      0.1414703654972616  -1.4703375146553483     0.3134894619564954
# rs75675305      chr6;109324353;T;C      50      13      1893    537     0.13105806668347084     0.532414709954525       0.7074774689897501  0.37524607180730757     0.7074774689897501
# rs7761290       chr6;109326621;T;G      40      8       2979    267     -0.6817584084935364     0.3095160694834158      0.19383573682018126 -1.2993153261065467     0.3134894619564954
result_table=args.result_table
result_table=pd.read_table(result_table,sep="\t",header=0)

### remove all the event/barcodes unrealted to our event of interest
index_event=index_event[index_event["mutation"].isin(result_table["event"])]

### get multiple vis info files and concatanate them.
input_list= args.input_vis_file.split(",")

vis_info = []
for filename in input_list:
    vis_info.append(pd.read_table(filename,sep="\t",names=["index","fragmentname","readname","chrid","start","end"]))
vis_info = pd.concat(vis_info)

### take only elements related to our events 
index_event=index_event.merge(vis_info,left_on="start_aln_UMI",right_on="index",how="inner")

chromsizes=pybed.chromsizes(os.path.basename(args.genome).split('.')[0])
#chromsizes=pybed.chromsizes("hg19")

### get rid of vis_info as we dont need it anymore.
del vis_info

for event in result_table["event"].to_list():
    snpid=list(result_table[result_table["event"] == event]["SNP"])[0] #NaN if its NA
    print(f'creating bedgraph for {event} {snpid}',file=sys.stderr)
    generate_bedgraph(index_event, 1, snpid, chromsizes, event,output_dir)
    print("WT done!",file=sys.stderr)
    generate_bedgraph(index_event, 2, snpid, chromsizes, event,output_dir)
    print("VAR done!",file=sys.stderr)

# # here for debug
# event = "chr1;153896370;T;C"
# snpid=list(result_table[result_table["event"] == event]["SNP"])[0] #NaN if its NA
# generate_bedgraph(index_event, 1, snpid, chromsizes, event,output_dir)


# event_type = 1
# reads=index_event[(index_event["mutation"] == event) & (index_event["type"] == event_type )][["chrid","start","end"]]
# reads=pybed.BedTool.from_dataframe(reads)
# reads=reads.genomecov(bg=True,g=chromsizes,scale=1/int(result_table[result_table["event"] == event]["ref_count"]))
# reads = reads.to_dataframe()