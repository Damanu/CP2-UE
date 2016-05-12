#usr/bin/python!
#Calculates the fouriertransformed of a dataset in a text file 
#and writes it into a txt file named fft_out.txt

import sys
import numpy as np
if len(sys.argv)!=4:
	print "No Argument given\nArguments: Filename xcolumn ycolumn"
	sys.exit()
FILE=sys.argv[1]
XCOL=int(sys.argv[2])
YCOL=int(sys.argv[3])

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
	ft=np.fft.fft(y)
	delx=x[1]-x[0]#calculate delta x
	omega_2=np.pi/delx#calculate nyquist frequency
	w=np.linspace(-omega_2,omega_2,len(ft))
	f=open("fft_out.txt","w")
	i=0
	for val in ft:
		f.write(str(w[i]))#print frequency
		f.write("\t")
		f.write(str(abs(val)))#print absolute value
		f.write("\t")
		f.write(str(np.real(val)))#print real value
		f.write("\t")
		f.write(str(np.imag(val)))#print imaginary value
		f.write("\n")
		i+=1

if __name__=="__main__":
	main()
