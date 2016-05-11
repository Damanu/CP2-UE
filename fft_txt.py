#usr/bin/python!

import sys
import numpy as np

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
	y=np.fft.fft(x)
	print(y)
	f=open("fft_out.txt","w")
	i=0
	for val in y:
		f.write(str(i))
		f.write("\t")
		f.write(str(abs(val)))
		f.write("\n")
		i+=1

if __name__=="__main__":
	main()
