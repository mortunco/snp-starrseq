import argparse


parser = argparse.ArgumentParser(description='usage -f css.FASTA -s mrsfast.SAM -o outputdir')
parser.add_argument('-f',required=True, help='fasta file input')
parser.add_argument('-s',required=True, help='sam file input')
parser.add_argument('-o',required=True, help='output dir')
parser.add_argument('-a',required=True, help='adapter length')

args = parser.parse_args()


#mrsfast_sam="merge.css.noheader.sam"
mrsfast_sam= args.s
total_fasta=args.f
output_dir=args.o

read_dict={}
clip_length=int(args.a)
with open(mrsfast_sam) as sam_in:
    while True:
        line=sam_in.readline().rstrip().split("\t")
        if line == [""]:
            break
        if line[0][0]=="@":
            continue
        #print(line)
        if line[2] not in read_dict:
            read_dict[line[2]]={"count" : 1, "loc" : [], "dr" :""}
            read_dict[line[2]]["count"]=1
            read_dict[line[2]]["loc"].append(int(line[3]))
            read_dict[line[2]]["dr"]+=line[0] 

        else:
            read_dict[line[2]]["count"] +=1
            read_dict[line[2]]["loc"].append(int(line[3]))
            read_dict[line[2]]["dr"] += line[0]
print("reading mrfast output done")


# In[4]:


#total_fasta="css.fasta"
seq_dict={}
with open(total_fasta) as fasta_in:
    while True:
        line=fasta_in.readline()
        if line == "":
            break
        readname=line.rstrip().split()[0][1::]
        seq=fasta_in.readline().rstrip()
        seq_dict[readname]=seq
print("readin css fasta done")


# In[5]:


#out_single_fasta=open(args.o+"/"+"single.fasta","w")
out_with_fasta=open(args.o+"/"+"withadapter.fasta","w")
out_without_fasta=open(args.o+"/"+"withoutadapter.fasta","w")
#out_multi_fasta=open(args.o+"/"+"multiadapter.fasta","w")

#bedtools getfasta -fi sample.fa -bed x3.bed
#with min-1 max+15 
#without min+15 max-1

for i,y in read_dict.items():
    #print(i,y)
    if y["count"]==1:
        pass
        #out_single_fasta.write(">" +i + "\n")
        #out_single_fasta.write(seq_dict[i])
    if y["count"]==2:
#         start=min(y["loc"])
#         end=max(y["loc"])
        out_with_fasta.write(">" + i + "\n")
        out_with_fasta.write(seq_dict[i][min(y["loc"])-1:max(y["loc"])+clip_length] +"\n")
        
        out_without_fasta.write(">" +i + "\n")
        out_without_fasta.write(seq_dict[i][min(y["loc"])-1+clip_length:max(y["loc"])-1] +"\n")
    if y["count"]>2:
        pass
        #out_multi_fasta.write(">" +i + "\n")
        #out_multi_fasta.write(seq_dict[i])

#out_single_fasta.close()
out_with_fasta.close()
out_without_fasta.close()
#out_multi_fasta.close()
print("writing files done")    

    


# In[ ]:




