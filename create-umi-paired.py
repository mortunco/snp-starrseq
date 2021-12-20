import sys

while True:
	mate1=sys.stdin.readline().rstrip().split("\t")
	#print(mate1)
	if mate1[0] == "":
		break
	mate2=sys.stdin.readline().rstrip().split("\t")
	#print(mate1)
	#print(mate2)

	if mate1[3] > mate2[3]:
		#first attempt: index, read barcode, readname
		#print("{}:{}-{}{}\t{}{}\t{}".format(mate1[2],mate2[3],mate2[9][0:3],mate1[9][-3::],mate2[9][0:12],mate1[9][-12::],mate2[0]))
		#pretty: index, read barcode, readname
		#print(f'{mate1[2]}:{mate2[3]}-{mate2[9][0:3]}{mate1[9][-3::]}\t{mate2[9][0:12]}{mate1[9][-12::]}\t{mate2[0]}')
		#index, read barcode, readname, chr, start, end 
		print(f'{mate1[2]}:{mate2[3]}-{mate2[9][0:3]}{mate1[9][-3::]}\t{mate2[9][0:12]}{mate1[9][-12::]}\t{mate2[0]}\t{mate2[2]}\t{int(mate2[3])-1}\t{int(mate2[3])+abs(int(mate2[8]))}')
	else:
		#print("{}:{}-{}{}\t{}{}\t{}".format(mate1[2],mate1[3],mate1[9][0:3],mate2[9][-3::],mate1[9][0:12],mate2[9][-12::],mate1[0]))
		print(f'{mate1[2]}:{mate1[3]}-{mate1[9][0:3]}{mate2[9][-3::]}\t{mate1[9][0:12]}{mate2[9][-12::]}\t{mate1[0]}\t{mate1[2]}\t{int(mate1[3])-1}\t{int(mate1[3])+abs(int(mate1[8]))}')
