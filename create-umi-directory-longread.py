import sys

while True:
	mate1=sys.stdin.readline().rstrip().split("\t")
	#print(mate1)
	if mate1[0] == "":
		break
	print("{}:{}-{}{}\t{}{}\t{}".format(mate1[2],mate1[3],mate1[9][0:3],mate1[9][-3::],mate1[9][0:12],mate1[9][-12::],mate1[0]))