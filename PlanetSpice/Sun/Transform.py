import numpy as np
import spiceypy as sp
from ..utc2et import utc2et
import DateTimeTools as TT
import os
from .. import Globals

lsk_path = Globals.SpicePath + '/lsk/naif0010.tls'
spk_kernel = Globals.SpicePath + '/bodies/de432s.bsp'
pck_kernel = Globals.SpicePath + '/bodies/pck00010.tpc'
hci_kernel = Globals.SpicePath + '/frames/tk/sunframes.tk'

def HCItoIAU_SUN(Date,ut,xi,yi,zi):
	'''
	Convert from HCI to IAU_SUN coordinates
	'''
	
	#create some arrays
	n = np.size(xi)

	if np.size(ut) == 1:
		ut = np.zeros(n,dtype='float32') + ut
	if np.size(xi) == 1:
		xi = np.array([xi])
		yi = np.array([yi])
		zi = np.array([zi])
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')

	#load kernels
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(pck_kernel)
	sp.furnsh(hci_kernel)


	#do each unique date to speed things up
	if np.size(Date) == 1 :
		et[0] = utc2et(Date,0.0)
		et = et[0] + ut*3600.0
	else:
		ud = np.unique(Date)
		for i in range(0,ud.size):
			use = np.where(Date == ud[i])[0]
			tmp = utc2et(ud[i],0.0)
			et[use] = tmp + ut[use]*3600.0

	#transform coords
	for i in range(0,n):
		rot = sp.pxform('HCI','IAU_SUN',et[i])
		x[i],y[i],z[i] = np.sum(rot*np.array([xi[i],yi[i],zi[i]]),axis=1)
		
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(pck_kernel)
	sp.unload(hci_kernel)

	return (x,y,z)



def HAEtoHCI(Date,ut,xi,yi,zi):
	'''
	Convert HAE to HCI coordinates
	
	'''
	
	
	
	#create some arrays
	n = np.size(xi)

	if np.size(ut) == 1:
		ut = np.zeros(n,dtype='float32') + ut
	if np.size(xi) == 1:
		xi = np.array([xi])
		yi = np.array([yi])
		zi = np.array([zi])
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')

	#load kernels
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(pck_kernel)
	sp.furnsh(hci_kernel)
	
	#do each unique date to speed things up
	if np.size(Date) == 1 :
		et[0] = utc2et(Date,0.0)
		et = et[0] + ut*3600.0
	else:
		ud = np.unique(Date)
		for i in range(0,ud.size):
			use = np.where(Date == ud[i])[0]
			tmp = utc2et(ud[i],0.0)
			et[use] = tmp + ut[use]*3600.0

	
	for i in range(0,n):
		rot = sp.pxform('ECLIPDATE','HCI',et[i])
		x[i],y[i],z[i] = np.sum(rot*np.array([xi[i],yi[i],zi[i]]),axis=1)
		
	#free kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(pck_kernel)
	sp.unload(hci_kernel)

	return (x,y,z)

