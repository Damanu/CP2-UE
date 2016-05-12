#usr/bin/python!
#Calculates the fouriertransformed of a dataset in a text file 
#and writes it into a txt file named fft_out.txt

import sys
import numpy as np
if len(sys.argv)!=3:
	print "No Argument given\nArguments: Filename xcolumn"
	sys.exit()
FILE=sys.argv[1]
XCOL=int(sys.argv[2])

def main():
	f=open(FILE,'r')
	x=[]
	for line in f:
		x.append(line.split()[XCOL])
	x=map(float,x)	
	f.close()
	print np.mean(x), " +- ", np.std(x)
	 

if __name__=="__main__":
	main()
