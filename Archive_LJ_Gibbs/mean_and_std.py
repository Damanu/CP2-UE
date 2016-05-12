#/usr/bin/python!

#Program for calculation of meanvalue and standard deviation of a dataset
#output in text file:mean	std

#libraries
import os
import sys
import numpy as np

if len(sys.argv) != 3:
	print "You need at least 2 Arguments:\n	1. Filename\n	2. Column (starting at 0)"
	sys.exit()
FILE=sys.argv[1]
COL=int(sys.argv[2])

def Read():
	try:
		data=[]
		x=[]
		text= open(FILE,'r')
		data=text.readlines()
		colname=data[0].split()[COL]
		for line in data[1:]:
			f=line.split()
			x.append(f[COL])
		x=map(float,x)
		text.close()
		return x
	except:
		print "something went wrong in the Read() subroutine"
		sys.exit()
	
def main():
	data=Read()	
	print np.mean(data)
	print np.std(data)
if __name__=="__main__":
	main()

