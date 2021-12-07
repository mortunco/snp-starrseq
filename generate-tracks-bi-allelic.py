import pysam


import argparse
import os

parser = argparse.ArgumentParser(description="Gets barcode of each snp event then separates two bam files (REF/ALT)")
parser.add_argument("-bdir","--barcode-dir",type=str,required=True,help="path to directory that contains barcode infor for each snp in different file")
parser.add_argument("-resf","--result-file",type=str,required=True,help="Yi result file")
parser.add_argument("-bread","--barcode-read",type=str,required=True,help="path to file that contains read and barcode info.")
parser.add_argument("-bam","--bam",type=str,required=True,help="short read hiseq raw mapping file")
parser.add_argument("-odir","--output-dir",type=str,default="subsampled_bams/",help="directory to dump ref and alt bam files.")


args = parser.parse_args()


target_dir=args.barcode_dir

barcode_dic={}

for i in os.listdir(target_dir):
	rsid=i.split("_reads.txt")[0]
	temp_rsid_barcodes=[]
	temp_barcode={}
	with open(target_dir + "/" + i) as barcode_type_file:
		for line in barcode_type_file:
			temp=line.rstrip().split("\t")
			temp_barcode[temp[0]]=temp[2]
	
	barcode_dic[rsid]=temp_barcode
	 


result_file = args.result_file


rsid_info={}
with open(result_file) as r_file:
	for line in r_file:
		temp=line.rstrip().split("\t")
		if temp[0] == "SNP":
			continue
		rsid_info[temp[0]]=[temp[1].split(";")[0], temp[1].split(";")[1]]
		 


barcode_read = args.barcode_read

read_barcode_dic={}

with open(barcode_read) as barcode_read_f:
	for line in barcode_read_f:
		temp=line.rstrip().split("\t")
		read_barcode_dic[temp[0]]=temp[1]


if os.path.exists(args.output_dir) == True:
	pass
else:
	os.mkdir(args.output_dir)
		
bamfile = pysam.AlignmentFile(args.bam, "rb")
for rsid in list(rsid_info.keys()):
	if rsid == "SNP":
		continue
	print(rsid)
	ref=[]
	alt=[]
	ref_bam = pysam.AlignmentFile(args.output_dir+ "/" + rsid_info[rsid][0] + "_" + rsid_info[rsid][1] + "_" + rsid + ".REF.bam" ,"wb", template=bamfile)
	
	alt_bam = pysam.AlignmentFile(args.output_dir+ "/" + rsid_info[rsid][0] + "_" + rsid_info[rsid][1] + "_" + rsid + ".ALT.bam" ,"wb", template=bamfile)
	
	for read in bamfile.fetch(rsid_info[rsid][0], int(rsid_info[rsid][1])-1000, int(rsid_info[rsid][1])+1000):
		if read.is_supplementary == True or read.is_unmapped == True:
				continue
		if (read_barcode_dic[read.qname] in barcode_dic[rsid]):
			if barcode_dic[rsid][read_barcode_dic[read.qname]] == "1":
				ref_bam.write(read)
			elif barcode_dic[rsid][read_barcode_dic[read.qname]] == "2":
				alt_bam.write(read)
			else:
				print("oh shit")
				print(rsid,read)
				
				
#pysam.index(args.output_dir+ "/" + rsid_info[rsid][0] + "_" + rsid_info[rsid][1] + "_" + rsid + ".REF.bam")
#pysam.index(args.output_dir+ "/" + rsid_info[rsid][0] + "_" + rsid_info[rsid][1] + "_" + rsid + ".ALT.bam")
bamfile.close()
ref_bam.close()
alt_bam.close()




