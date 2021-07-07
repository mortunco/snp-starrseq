import pysam
import warnings

### format is going to be as follows
### OUTPUT FORMAT IS 1 BASED.
cigar_identifiers=["M","I","D","N","S",'H',"P","=","X","B"]


import argparse

# required arg

parser = argparse.ArgumentParser()
   
parser.add_argument('--bam_input', required=True)
args = parser.parse_args()



samfile = pysam.AlignmentFile(args.bam_input, "rb")
for read in samfile.fetch():
	read_cigar=[]
	if read.is_supplementary or read.is_unmapped:
		continue
	if read.get_tag("NM") == 0:
		print("\t".join([read.seq[0:12] + read.seq[-12::],read.reference_name + ":" + str(read.reference_start+1) + "-" + read.seq[0:3] + read.seq[-3::],read.qname,"WT","NA" ,"NA","NA","NA","PASS"]))
		continue
	for cigar_pair in read.cigartuples:
		read_cigar.append(cigar_identifiers[cigar_pair[0]] * cigar_pair[1])
	ins_check=False
	del_check=False
	# print(read.qname)
	# print(read.seq)
	# print("xddd")
	del_corrector=0
	for i,(ref,cigar) in enumerate(zip(read.get_aligned_pairs(with_seq=True),"".join(read_cigar))):
		if cigar == "D":
			#ref is something like this
			#(0, s15905, 'T')
			temp=["None",ref[0],ref[1],ref[2],cigar]
			del_corrector +=1
		else:
			temp=[read.seq[i-del_corrector],ref[0],ref[1],ref[2],cigar]
		#print(ref)
		#print(temp)
		### temp ###
		# [READ BASE, 
		# POSITION OF THAT BASE ON READ, 
		# POSITION OF REFERENCE GENOME, 
		# REFERENCE GENOME BASE, 
		# POSITION CIGAR]
		if cigar == "M":
			### insertion calling ####
			if ins_check == False and del_check == False:
				previous_temp=temp

			if ins_check == True:
				ins_check = False
				#variantline=[read.qname,"ins",str(temp[2]-len(insertion) + 1),insertion[0],insertion]
				variantline=[read.seq[0:12] + read.seq[-12::],read.reference_name + ":" + str(read.reference_start+1) + "-" + read.seq[0:3] + read.seq[-3::],read.qname,"ins",read.reference_name, str(previous_temp[2]+1),previous_temp[3].upper(),previous_temp[3].upper()+insertion.upper()]
				if "H" in read_cigar:
					print("\t".join(variantline+["SupMappingBarcode"]) )
				else:
					print("\t".join(variantline+["PASS"]))

			### deletion calling ###
			if del_check == True:
				del_check=False
				variantline=[read.seq[0:12] + read.seq[-12::],read.reference_name + ":" + str(read.reference_start+1) + "-" + read.seq[0:3] + read.seq[-3::],read.qname,"del",read.reference_name,str(previous_temp[2]+1),previous_temp[3].upper()+deletion.upper(),previous_temp[3].upper()]
				if "H" in read_cigar:
					print("\t".join(variantline+["SupMappingBarcode"]) )
				else:
					print("\t".join(variantline+["PASS"]) ) 
			### snp calling ###
			if temp[0] == temp[3]:
				pass
			else:
				variantline=[read.seq[0:12] + read.seq[-12::],read.reference_name + ":" + str(read.reference_start+1) + "-" + read.seq[0:3] + read.seq[-3::],read.qname,"snp",read.reference_name,str(temp[2] + 1),temp[3].upper(),temp[0]]	
				if "H" in read_cigar:
					print("\t".join(variantline+["SupMappingBarcode"]) )
				else:
					print("\t".join(variantline + ["PASS"] ) )
		elif cigar == "I":
			if ins_check == False:
				### last match sequence ###
				#insertion=read.seq[i-1] + temp[0]
				insertion= temp[0]
				ins_check=True
			else:
				insertion=insertion+temp[0]
		elif cigar == "D":
			if del_check == False:
				#deletion = read.seq[i-1] + temp[3]
				deletion =  temp[3]
				del_check=True
			else:
				deletion = deletion + temp[3]
		elif cigar == "S":
			pass
		elif cigar == "H":
			pass
		else:
			warnings.warn("Off cound cigar type found --> {}. Becareful with the cigar {}".format(read_cigar,cigar))



	


#print(len(read.get_aligned_pairs(with_seq=True)), len(read.seq),len("".join(read_cigar)))
