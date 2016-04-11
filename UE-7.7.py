#/usr/bin/python!

import numpy as np
import sys
import matplotlib.pyplot as plt
#-----------Help----------------
if len(sys.argv)==1:
	print "Arguments: input file, dr, N"
	sys.exit()

#-----------Globals------------
PI=np.pi
PATH="./"
N=int(sys.argv[1])
INFILE="data.txt"
roh=0.85

#----------------Routines-----------------
def Read():
	data=[]
	try:
		datax=[]
		datay=[]
		text= open(PATH+INFILE,'r')
		data=text.readlines()
		print data[1:]
		for x in data[1:]:
			x_=x.split()
			print x_
			datax.append(x_[0])	
			datay.append(x_[1])	
		text.close()
		print "in"
		return datax,datay
	except: 
		print 'something went wrong at subprogram Read() while reading data file'
		sys.exit()

def integrand_g(g,r,k,roh):
	f=4*PI*roh/k*r*np.sin(k*r)*(g-1)
	return f

def integrand_s(S,r,k,roh):
	f=1./(2.*PI**2*roh*r)*k*np.sin(k*r)*(S-1)
	return f

#--------------------Main---------------------------------
def main():
	S=[]
	g=[]

	dr=float(sys.argv[2])
	dk=20./N		#k is in units of 1/sig 
	r=np.linspace(0,(N-1)*dr,N)
#	print r
	[x,y]=Read()
	y=np.array(map(float,y))
	x=np.array(map(float,x))
	
	k=dk
#	print integrand_g(y,r,k,roh)
	while k < N*dk:
		S.append(1+np.trapz(integrand_g(y,r,k,roh),x=None, dx=dr))
		k+=dk
#		print k
	S=np.array(S)	
	r=dr
	while r < N*dr:
		g.append(1+np.trapz(integrand_s(S,r,k,roh),x=None, dx=dk))
		r+=dr
#		print r
	
#	print "S(k)=",S
#	print "g(r)=",g
#	plt.plot(np.linspace(0,N*dk,len(S)),S)
	print S
	plt.plot(x,y)
	plt.plot(np.linspace(0,N*dk,len(S)),S)
	plt.show()
if __name__=="__main__":
	main()
