import argparse
import string
parser = argparse.ArgumentParser(
    description="gets collapsed reads and matches long reads of the two independent runs")
parser.add_argument("--forward_fastq",
  type=str,
  required=True,
  help="Input forward fastq")
parser.add_argument("--reverse_fastq",
  type=str,
  required=True,
  help="Input reverse fastq")

parser.add_argument("--prefix",
  type=str,
  required=True,
  help="run prefix")
args = parser.parse_args()


def rev_compl(st):
    ### from https://stackoverflow.com/a/60070333 ###
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "N" :"N"}
    return "".join(nn[n] for n in reversed(st))

forward=open(args.forward_fastq,"r")
reverse=open(args.reverse_fastq,"r")

barcode_to_cid = {} 
cid_to_readinfo= {}

read_counter=0
while True:
  forward_readname = forward.readline().strip()
  forward_seq = forward.readline().strip()
  forward_qual = forward.readline().strip()
  forward_qual = forward.readline().strip()
  read_counter+=1
  if not forward_readname:
   break
  reverse_readname = reverse.readline().strip()
  reverse_seq = reverse.readline().strip()
  reverse_qual = reverse.readline()
  reverse_qual = reverse.readline().strip()
  read_counter+=1
  barcode= forward_seq[0:12] + rev_compl(reverse_seq[0:12])
  #print(barcode)
  cid=forward_readname.lstrip("@").rstrip().split("\t")[0] ### Forward term is related to merged files. Not longshort or short long. ALSO, split should be space if there merged fastqs are piped in to seqtk
  #print([cid,barcode])
  if barcode not in barcode_to_cid:
    barcode_to_cid[barcode]=[cid]
  else:
    barcode_to_cid[barcode].append(cid)
  cid_to_readinfo[cid]={"forward" : {"seq": forward_seq, "qual": forward_qual}, "reverse": {"seq": reverse_seq, "qual":  reverse_qual}}
    #cid_to_readinfo[cid]=[forward_seq, forward_qual, reverse_seq, reverse_qual]

forward.close()
reverse.close()
#print(len(barcode_to_cid["GGAGTTGCAAGTTATGTCCCCAAA"]))
#print(cid_to_readinfo["12longshort"])
#print(cid_to_readinfo["39shortlong"])

orphan=open(args.prefix+"/orphan.fastq","w")
clustered=open(args.prefix+"/clustered.interleaved.fastq","w")
problematic_samesame=open(args.prefix+"/problematic_samesame.interleaved.fastq","w")
problematic_multiple=open(args.prefix+"/problematic_multiple.interleaved.fastq","w")
master_file=open(args.prefix+"/master-barcode-cid.txt","w")

#print(barcode_to_cid)
for barcode,clusters in barcode_to_cid.items():
  if len(clusters) == 1:
    orphan.write("@"+barcode+"\n")
    orphan.write(cid_to_readinfo[clusters[0]]["forward"]["seq"] +"\n")
    orphan.write("+\n")
    orphan.write(cid_to_readinfo[clusters[0]]["forward"]["qual"] +"\n")
    
    orphan.write("@"+barcode+"\n")
    orphan.write(cid_to_readinfo[clusters[0]]["reverse"]["seq"] +"\n")
    orphan.write("+\n")
    orphan.write(cid_to_readinfo[clusters[0]]["reverse"]["qual"] +"\n")
    master_file.write(barcode+"\t"+clusters[0]+ "\t" "orphan"+"\n")
  elif len(clusters) == 2:
    #print(clusters)
    check_names_same=[]
    for i in clusters:
      if "longshort" in i:
        check_names_same.append(True)
      else:
        check_names_same.append(False)
    if (check_names_same == [False,False]) or (check_names_same == [True,True]):
      for index,cid in enumerate(clusters):
        master_file.write(barcode+"\t"+cid+"\t"+"problematic_same_same"+"\n")  
        problematic_samesame.write("@"+barcode+".problematic_"+str(index)+"\n")
        problematic_samesame.write(cid_to_readinfo[cid]["forward"]["seq"] +"\n")
        problematic_samesame.write("+\n")
        problematic_samesame.write(cid_to_readinfo[cid]["forward"]["qual"] +"\n")
        problematic_samesame.write("@"+barcode+".problematic_"+str(index)+"\n")
        problematic_samesame.write(cid_to_readinfo[cid]["reverse"]["seq"] +"\n")
        problematic_samesame.write("+\n")
        problematic_samesame.write(cid_to_readinfo[cid]["reverse"]["qual"] +"\n")
    else:
      master_file.write(barcode+"\t"+clusters[0]+"\t"+"clustered"+"\n")
      master_file.write(barcode+"\t"+clusters[1]+"\t"+"clustered"+"\n")
      clustered.write("@"+barcode+"\n")
      clustered.write(max(cid_to_readinfo[clusters[0]]["forward"]["seq"],cid_to_readinfo[clusters[1]]["forward"]["seq"],key=len) + "\n")
      clustered.write("+\n")
      clustered.write(max(cid_to_readinfo[clusters[0]]["forward"]["qual"],cid_to_readinfo[clusters[1]]["forward"]["qual"],key=len) + "\n")
      clustered.write("@"+barcode+"\n")
      clustered.write(max(cid_to_readinfo[clusters[0]]["reverse"]["seq"],cid_to_readinfo[clusters[1]]["reverse"]["seq"],key=len) + "\n")
      clustered.write("+\n")
      clustered.write(max(cid_to_readinfo[clusters[0]]["reverse"]["qual"],cid_to_readinfo[clusters[1]]["reverse"]["qual"],key=len) + "\n")
  else:
    #print(clusters)
    #print(len(clusters))
    for i in clusters:
      #print(i)
      master_file.write(barcode+"\t"+i+"\t"+"problematic_multiple"+"\n")  
      problematic_multiple.write("@"+barcode+"_"+i+"\n")
      problematic_multiple.write(cid_to_readinfo[i]["forward"]["seq"] +"\n")
      problematic_multiple.write("+\n")
      problematic_multiple.write(cid_to_readinfo[i]["forward"]["qual"] +"\n")
      problematic_multiple.write("@"+barcode+"_"+i+"\n")
      problematic_multiple.write(cid_to_readinfo[i]["reverse"]["seq"] +"\n")
      problematic_multiple.write("+\n")
      problematic_multiple.write(cid_to_readinfo[i]["reverse"]["qual"] +"\n")

      #print(cid_to_readinfo[i]["forward"]["seq"],cid_to_readinfo[i]["reverse"]["seq"])

    #raise ValueError('more than two clusters with same barcode ({}) danger danger danger'.format(barcode))
    
orphan.close()
reverse.close()
problematic_samesame.close()
problematic_multiple.close()


# for barcode,clusters in barcode_to_cid.items():
#   for i in clusters:
#     master_file.write(barcode+"\t"+i+"\n")

master_file.close()
#print(read_counter)

