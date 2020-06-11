import numpy as np
import spiceypy as sp
from ..utc2et import utc2et
from scipy.interpolate import InterpolatedUnivariateSpline
import os
import DateTimeTools as TT
from ..Tools.FileSearch import FileSearch

lsk_path = Globals.SpicePath + '/lsk/naif0010.tls'
spk_kernel = Globals.SpicePath + '/bodies/de432s.bsp'
pck_kernel = Globals.SpicePath + '/bodies/pck00010.tpc'

hci_kernel = Globals.SpicePath + '/frames/tk/sunframes.tk'
vso_kernel = Globals.SpicePath + '/VEX/tk/VenusVSO.tk'

def ListVenusSPK(Date):
	'''
	Lists the kernels for VEX
	
	'''
	if np.size(Date) == 1:
		yymm = np.array([(Date%1000000)/100])
	else:
		yymm = np.array((Date%1000000)/100)
	
	yymm = np.unique(yymm)
	
	yy = yymm/100
	mm = yymm%100
	
	#pad with an extra month either side
	if mm[0] > 1:
		mm = np.append(mm[0]-1,mm)
		yy = np.append(yy[0],yy)
	else:
		mm = np.append(12,mm)
		yy = np.append(yy[0]-1,yy)
	
	if mm[-1] < 12:
		mm = np.append(mm,mm[-1]+1)
		yy = np.append(yy,yy[-1])
	else:
		mm = np.append(mm,1)
		yy = np.append(yy,yy[-1]+1)
	yymm = yy*100 + mm
	
	#check for months which don't exist
	use = np.where((yymm >= 604) & (yymm <= 1412))[0]
	yymm = yymm[use].astype('int64')
	yymmdd = (yymm*100 + 1)*1000000
	if yymmdd[0] == 60401000000:
		yymmdd[0] = 60409211524
		
	#now search for files
	n = np.size(yymm)
	path =  Globals.SpicePath + '/VEX/spk/'
	files = np.zeros(n,dtype='S59')
	for i in range(0,n):
		tmp = FileSearch(path,'*{:012d}*.BSP'.format(yymmdd[i]))
		files[i] = path + tmp[0]
	
	return files


		
def PosVSO(Date,ut):
	'''
	VEX position in VSO coordinates
	
	'''
	
	#create the output arrays
	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
		
	#load the relevant kernels
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	VEXspk = ListVenusSPK(date)
	for vk in VEXspk:
		sp.furnsh(vk)
	sp.furnsh(pck_kernel)
	sp.furnsh(vso_kernel)	
	
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
			
	#get the positions for each date/time
	pos,lt = sp.spkpos('VEX',et,'VENUSVSO','NONE','VENUS')
	x = pos.T[0]
	y = pos.T[1]
	z = pos.T[2]		
	
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	for vk in VEXspk:
		sp.unload(vk)
	sp.unload(pck_kernel)
	sp.unload(vso_kernel)

	return (x,y,z)


def PosHCI(Date,ut):
	'''
	VEX position in HCI coordinates
	
	'''
	
	#create the output arrays
	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
		
	#load the relevant kernels
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	VEXspk = ListVenusSPK(date)
	for vk in VEXspk:
		sp.furnsh(vk)
	sp.furnsh(pck_kernel)
	sp.furnsh(vso_kernel)	
	
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
			
	#get the positions for each date/time
	pos,lt = sp.spkpos('VEX',et,'HCI','NONE','SUN')
	x = pos.T[0]
	y = pos.T[1]
	z = pos.T[2]		
	
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	for vk in VEXspk:
		sp.unload(vk)
	sp.unload(pck_kernel)
	sp.unload(vso_kernel)

	return (x,y,z)






def CarringtonLongitude(Date,ut):
	
	#create output array
	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	lon = np.zeros(n,dtype='float64')
	if n == 1:
		ut = np.array([ut])
		
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
	
	#get the longitudes
	pos,lt = sp.spkpos('VEX',et,'IAU_SUN','NONE','SUN')
	lon = np.arctan2(pos.T[1],pos.T[0])
	
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(pck_kernel)	
	sp.unload(hci_kernel)	

	return lon
