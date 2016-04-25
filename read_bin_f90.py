#/usr/bin/python!

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.io import FortranFile
BinaryData=sys.argv[1]
def read():
	f=FortranFile(BinaryData,'r')
	n_vals=f.read_record('i4,i4,i4,i4,i4,(12,)i4')
	d_vals=f.read_reals('f4')
	xyz_vals=f.read_record('(500,)f4,(500,)f4,(500,)f4')
	a_vals=f.read_reals('i8')
	ndr_vals=f.read_ints('i4')
	ag_vals=f.read_record('(3,)i8')
	f.close()
	return n_vals,d_vals,xyz_vals,a_vals,ndr_vals,ag_vals



def main():
	(n,d,xyz,a,ndr,ag)=read()
#	for i in ag:
#		print i[0],"\t",i[1],"\t",i[2]
	i=0
	while i<len(xyz[0]):
		print xyz[0][i],"\t",xyz[1][i],"\t",xyz[2][i]
		i+=1
#	plt.plot(np.array(0,1,len(f)),f)
#	plt.show()
if __name__=="__main__":
	main()
