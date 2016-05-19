#/usr/bin/python!

import numpy as np
import sys
import matplotlib.pyplot as plt
BinaryData=sys.argv[1]
bin_UE7_1="/home/emanuel/Git/Hub/CP2-UE/lj-canonical/mclj_out.dat"
dr=float(sys.argv[2])
#---------------datatypes-------------------
n=500

rho=0.84
#-------------------------------------------

#read data-----------------------------------------
def read():
	n=np.fromfile(BinaryData,dtype=np.dtype("i4"),count=1)[0]
	dr=np.fromfile(BinaryData,dtype=np.dtype("(8)i4,f8"),count=1)["f1"]
	rho=np.fromfile(BinaryData,dtype=np.dtype("(8)i4,f8,f8,f8"),count=1)["f3"]
	ndr=int(0.5*(n/rho)**(1/3.)/dr)
	print ndr
	dt_UE8_1=np.dtype([("n","i4"),("ncor","i4"),("nt","i4"),("ntaskip","i4"),("ntcskip","i4"),("ntjob","i4"),("ntorig","i4"),("ntprint","i4"),("dr","f8"),("dt","f8"),("rho","f8"),("x","f8",n),("y","f8",n),("z","f8",n),("vx","f8",n),("vy","f8",n),("vz","f8",n),("ak","f8"),("ak2","f8"),("au","f8"),("aw","f8"),("ag","i8",ndr)])
	data=np.fromfile(BinaryData,dtype=dt_UE8_1,count=1)
	return data_,data,



def main():
	E=[]
	data=read()
	ndr=int(0.5*(data["n"]/data["rho"])**(1/3.)/data["dr"])
	print "Epot= ",data["au"]*data["ntaskip"]/data["ntjob"],"p=",rho*(2.*data["ak"]*data["ntaskip"]/data["ntjob"]-data["aw"]*data["ntaskip"]/data["ntjob"])/(3*data["n"])
	i=0
	for x in ag[0]:
		g.append(x*data["ntaskip"]/(data["ntjob"]*2*np.pi*(i*dr)**2*dr*rho*n))
		i+=1
	i=0
	print "Epot= ",data_["au"]/data_["ntjob"],"p= ", data_["rho"]*(data_["t"]-data_["aw"]/(3.0*data_["ntjob"]));
	ag_[0][0]=0
	for x in ag_[0][1:]:
		g_.append(x/(data_["ntjob"]*2*np.pi*(i*data_["dr"])**2*data_["dr"]*rho*n))
		i+=1
	plt.plot(np.linspace(0,float(ndr)*float(dr),len(g)),g)
	plt.plot(np.linspace(0,float(ndr)*float(dr),len(g_)),g_)
	plt.show()
if __name__=="__main__":
	main()
