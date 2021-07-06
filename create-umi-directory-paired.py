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
		print("{}:{}-{}{}\t{}{}\t{}".format(mate1[2],mate2[3],mate2[9][0:3],mate1[9][-3::],mate2[9][0:12],mate1[9][-12::],mate2[0]))
	else:
		print("{}:{}-{}{}\t{}{}\t{}".format(mate1[2],mate1[3],mate1[9][0:3],mate2[9][-3::],mate1[9][0:12],mate2[9][-12::],mate1[0]))
