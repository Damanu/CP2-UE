#usr/bin/python!
#Calculates the fouriertransformed of a dataset in a text file 
#and writes it into a txt file named fft_out.txt

import sys
import numpy as np
if len(sys.argv)!=4:
	print "No Argument given\nArguments: Filename xcolumn ycolumn"
	sys.exit()
FILE=sys.argv[1]
YCOL=int(sys.argv[3])
XCOL=int(sys.argv[2])

def main():
	f=open(FILE,'r')
	x=[]
	y=[]
	for line in f:
		y.append(line.split()[YCOL])
		x.append(line.split()[XCOL])
	x=map(float,x)	
	y=map(float,y)	
	f.close()
	delx=x[1]-x[0]#calculate delta x
	intf=1/3.*np.trapz(y,x=x)
	print intf
	f=open("int_out.txt","w")
	f.write(str(intf))#print integral

if __name__=="__main__":
	main()
