#/usr/bin/python!

import numpy as np
import sys
import matplotlib.pyplot as plt
#-----------Help----------------
if len(sys.argv)==1:
	print "Arguments: N,rho"
	sys.exit()

#-----------Globals------------
PI=np.pi
PATH="./"
N=int(sys.argv[1])
INFILE="data.txt"
rho=float(sys.argv[2])

#----------------Routines-----------------
def Read():
	data=[]
	try:
		x=[]
		y=[]
		text= open(PATH+INFILE,'r')
		data=text.readlines()
		for line in data[1:]:
			f=line.split()
			x.append(f[0])	
			y.append(f[1])
		x=map(float,x)	
		y=map(float,y)	
		text.close()
		return x,y
	except: 
		print 'something went wrong at subprogram Read() while reading data file'
		sys.exit()
#--------------------Main---------------------------------
def main():
	[r,g]=Read()
	imax=g.index(max(g))
	g_=g[imax:]
	imin2=g_.index(min(g_))+imax-1
	r=np.array(r[:imin2])
	g=np.array(g[:imin2])
	n_r=np.trapz(4*PI*rho*r*r*g,x=None, dx=r[1]-r[0])
	print n_r
if __name__=="__main__":
	main()
