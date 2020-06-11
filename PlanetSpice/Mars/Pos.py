import numpy as np
import spiceypy as sp
import DateTimeTools as TT
import os
from scipy.interpolate import InterpolatedUnivariateSpline
from ..utc2et import utc2et
from .. import Globals
from ..Tools.ListDates import ListDates
from ..Tools.ContUT import ContUT

lsk_path = Globals.SpicePath + '/lsk/naif0010.tls'
spk_kernel = Globals.SpicePath + '/bodies/de432s.bsp'
pck_kernel = Globals.SpicePath + '/bodies/pck00010.tpc'
hci_kernel = Globals.SpicePath + '/frames/tk/sunframes.tk'


def PosHCI(Date,ut):
	'''
	HCI position at set times
	
	'''	

	#get output arrays
	n = np.size(ut)
	if np.size(ut) == 1:
		ut = np.array([ut])
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
			
	#get the positions
	pos,lt = sp.spkpos('4',et,'HCI','NONE','SUN')
	x = pos.T[0]
	y = pos.T[1]
	z = pos.T[2]	
	
	#unload kernels		
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(pck_kernel)
	sp.unload(hci_kernel)

	return (x,y,z)

def PosHAE(Date,ut):
	'''
	HAE position at set times
	
	'''	

	#get output arrays
	n = np.size(ut)
	if np.size(ut) == 1:
		ut = np.array([ut])
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
	
	#get the positions
	pos,lt=sp.spkpos('4',et,'ECLIPDATE','NONE','SUN')
	x = pos.T[0]
	y = pos.T[1]
	z = pos.T[2]	
	
	#unload kernels		
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(pck_kernel)
	sp.unload(hci_kernel)

	return (x,y,z)



def SaveCarringtonRotations(Date0=19500101,Date1=20500101):
	'''
	Create a list of times when Mars is at a solar longitude of 0,
	defining the start of a Carrington rotation.
	
	'''

	#name the file to save the data in
	outpath = Globals.OutputPath + 'Mars/'
	if not os.path.isdir(outpath):
		os.system('mkdir -pv '+outpath)
	outfile = outpath + '0long.dat'
		
	#list the dates between date0 and date1
	dates = ListDates(Date0,Date1)
	nd = dates.size
	nt = nd*24
	
	#load the kernels
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(pck_kernel)
	sp.furnsh(hci_kernel)		
	
	#get the ets
	et0 = utc2et(dates[0],0.0)
	et = et0 + np.arange(nt,dtype='float64')*3600
	
	#get the array of longitudes
	lon = np.zeros(nt,dtype='float64')
	for i in range(0,nt):
		print('\rFinding longitude {0} of {1}'.format(i+1,nt),end='')
		pos,lt = sp.spkpos('4',et[i],'IAU_SUN','NONE','SUN')
		lon[i] = np.arctan2(pos[1],pos[0])
	print('')
		
			
	#find the bits where we cross 0
	n=0
	for i in range(0,nt-1):
		print('\rFinding 0s {:6.2f}%'.format(100*np.float32((i+1))/(nt-1)),end='')
		if lon[i] > 0 and lon[i+1] <= 0:
			if i == 0:
				inds = np.array([3,2,1,0])
			elif i >= nt-2:
				inds = np.array([nt-1,nt-2,nt-3,nt-4])
			else:
				inds = np.arange(4)[::-1]+i-1
			f = InterpolatedUnivariateSpline(lon[inds],et[inds])
			if n == 0:
				ets = np.array([f(0.0)],dtype='float64')
			else:
				ets = np.append(ets,f(0.0))
			n+=1

	print('')
	
	
	#store them in a file
	f=open(outfile,'w')
	for i in range(0,n):
		print('\rSaving time {0} of {1}'.format(i+1,n),end='')
		s = sp.timout(ets[i],"YYYYMMDD HRMNSC ::RND",32)
		date,ut = s.split()
		date = np.int32(date)
		ut = TT.HHMMtoDec(np.int32(ut),ss=True)
		utc = ContUT(np.array([date]),np.array([ut]))[0]
		outstr = '{:08d} {:f} {:f}\n'.format(date,ut,utc)
		f.write(outstr)
	f.close()
	print('')

	#unload kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(pck_kernel)	
	sp.unload(hci_kernel)	



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
	pos,lt = sp.spkpos('4',et,'IAU_SUN','NONE','SUN')
	lon = np.arctan2(pos.T[1],pos.T[0])
	
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(pck_kernel)	
	sp.unload(hci_kernel)	

	return lon



def ReadCarringtonRotations():
	
	path = Globals.OutputPath + 'Mars/'
	fname = path + '0long.dat'
	
	return pf.ReadASCIIData(fname,Header=False,dtype=dtypecarr)


