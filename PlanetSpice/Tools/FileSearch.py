import fnmatch as fnm
import os as os
import numpy as np

def FileSearch(dirname,fname):
	if os.path.isdir(dirname) == False:
		return np.array([])
		
	files=np.array(os.listdir(dirname))
	files.sort()
	matches=np.zeros(np.size(files),dtype='bool')
	for i in range(0,np.size(files)):
		if fnm.fnmatch(files[i],fname):
			matches[i]=True
	
	good=np.where(matches == True)[0]
	if np.size(good) == 0:
		return np.array([])
	else:
		return np.array(files[good])
