#/usr/bin/python!

import numpy as np
import sys
from scipy.io import FortranFile
BinaryData=sys.argv[1]
def read():
	filearray=FortranFile(BinaryData,'r')
	print filearray
	f=filearray.read_reals()
	filearray.close()
	return f



def main():
	f=read()
	print f
if __name__=="__main__":
	main()
