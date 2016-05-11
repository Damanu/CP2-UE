#usr/bin/python!

import sys

FILE=sys.argv[1]
COL=int(sys.argv[2])

def main():
	f=open(FILE,'r')
	x=[]
	for line in f:
		x.append(line.split()[COL])
	x=map(float,x)	
	print x
	f.close()
if __name__=="__main__":
	main()
