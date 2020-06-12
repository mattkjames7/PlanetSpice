import numpy as np
import spiceypy as sp
import DateTimeTools as TT
import os
from scipy.interpolate import InterpolatedUnivariateSpline
from ..utc2et import utc2et
from .. import Globals
import PyFileIO as pf
from ..Tools.ListDates import ListDates
from ..Tools.ContUT import ContUT
import PyFileIO as pf
import RecarrayTools as RT

lsk_path = Globals.SpicePath + '/lsk/naif0010.tls'
spk_kernel = Globals.SpicePath + '/bodies/de432s.bsp'
pck_kernel = Globals.SpicePath + '/bodies/pck00010.tpc'
hci_kernel = Globals.SpicePath + '/frames/tk/sunframes.tk'

#a common dtype used for storing position
dtype = [	('Date','int32'),
			('ut','float32'),
			('utc','float64'),
			('xHCI','float64'),
			('yHCI','float64'),
			('zHCI','float64'),
			('xIAU_SUN','float64'),
			('yIAU_SUN','float64'),
			('zIAU_SUN','float64'),
			('Rsun','float64'),
			('LatHCI','float32'),
			('LonHCI','float32'),
			('LatIAU_SUN','float32'),
			('LonIAU_SUN','float32')]



#carrington rot dtype
dtypecarr = [	('Date','int32'),
				('ut','float32'),
				('utc','float64')]
AU = 1.496e8

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
	pos,lt = sp.spkpos('VENUS',et,'HCI','NONE','SUN')
	pos = np.array(pos)
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
	pos,lt=sp.spkpos('VENUS',et,'ECLIPDATE','NONE','SUN')
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
	Create a list of times when Venus is at a solar longitude of 0,
	defining the start of a Carrington rotation.
	
	'''

	#name the file to save the data in
	outpath = Globals.OutputPath + 'Venus/'
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
		pos,lt = sp.spkpos('VENUS',et[i],'IAU_SUN','NONE','SUN')
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
	pos,lt = sp.spkpos('VENUS',et,'IAU_SUN','NONE','SUN')
	lon = np.arctan2(pos.T[1],pos.T[0])
	
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(pck_kernel)	
	sp.unload(hci_kernel)	

	return lon



def ReadCarringtonRotations():
	
	path = Globals.OutputPath + 'Venus/'
	fname = path + '0long.dat'
	
	return pf.ReadASCIIData(fname,Header=False,dtype=dtypecarr)




def PosIAU_SUN(Date,ut):
	'''
	Get Venus' position in IAU_SUN coordinates, wher Z is along the 
	Sun's rotational axis, X and Y rotate with the Sun.
	
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
	pos,lt = sp.spkpos('VENUS',et,'IAU_SUN','NONE','SUN')
	pos = np.array(pos)
	x = pos.T[0]
	y = pos.T[1]
	z = pos.T[2]	
	
	#unload kernels	
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(pck_kernel)
	sp.unload(hci_kernel)

	return (x,y,z)


def SavePos(Date0=19500101,Date1=20500101):
	'''
	Save Venus' position for every date between date0 and date1
	'''
	outpath = Globals.OutputPath + 'Venus/VenusPos/'
	if not os.path.isdir(outpath):
		os.system('mkdir -pv '+outpath)

	
	date = Date0
	ut = np.arange(24.0)
	while date <= Date1:
		print('Saving date {:d}'.format(date))
		fname = outpath + '{:08d}.bin'.format(date)
		
		#calculate the postions
		x,y,z = PosHCI(date,ut)
		x2,y2,z2 = PosIAU_SUN(date,ut)
		
		#create the output array
		data = np.recarray(24,dtype=dtype)
		
		#fill it
		data.Date = date
		data.ut = ut
		data.utc = ContUT(data.Date,data.ut)
		
		data.xHCI = x/AU
		data.yHCI = y/AU
		data.zHCI = z/AU
		data.xIAU_SUN = x2/AU
		data.yIAU_SUN = y2/AU
		data.zIAU_SUN = z2/AU
		
		#calculate some things
		data.Rsun = np.sqrt(x**2 + y**2 + z**2)/AU
		xyHCI = np.sqrt(data.xHCI**2 + data.yHCI**2)
		xyIAU_SUN = np.sqrt(data.xIAU_SUN**2 + data.yIAU_SUN**2)
		data.LatHCI = np.arctan2(data.zHCI,xyHCI)*180.0/np.pi
		data.LonHCI = np.arctan2(data.yHCI,data.xHCI)*180.0/np.pi
		data.LatIAU_SUN = np.arctan2(data.zIAU_SUN,xyIAU_SUN)*180.0/np.pi
		data.LonIAU_SUN = np.arctan2(data.yIAU_SUN,data.xIAU_SUN)*180.0/np.pi		
		
		#save file
		RT.SaveRecarray(data,fname)
		
		#update date
		date = TT.PlusDay(date)
		
		

def ReadPosDate(Date):
	outpath = Globals.OutputPath + 'Venus/VenusPos/'
	fname = outpath + '{:08d}.bin'.format(Date)
	if not os.path.isfile(fname):
		return np.recarray(0,dtype=dtype)
	return RT.ReadRecarray(fname,dtype)
	
def ReadPos(Date0,Date1):
	
	Dates = ListDates(Date0,Date1)
	nd = Dates.size
	
	n = 0
	path = Globals.OutputPath + 'Venus/VenusPos/'
	for i in range(0,nd):
		print('\rCounting records: {:d}'.format(n),end='')
		fname = path + '{:08d}.bin'.format(Dates[i])
		if os.path.isfile(fname):
			f = open(fname,'rb')
			n += np.fromfile(f,dtype='int32',count=1)[0]
			f.close()
	print('\rCounting records: {:d}'.format(n))

	out = np.recarray(n,dtype=dtype)
	
	p = 0
	for i in range(0,nd):
		print('\rReading Date {:08d} ({:05.1f}%)'.format(Dates[i],100.0*np.float(i+1)/nd), end=' ')
		tmp = ReadPosDate(Dates[i])
		out[p:p+tmp.size] = tmp
		p += tmp.size
	print('')
	
	return out



def CombinePos(Date0=19500101,Date1=20500101):
	data = ReadPos(Date0,Date1)
	fname = Globals.OutputPath + 'Venus/VenusPos.bin'
	RT.SaveRecarray(data,fname) 
	

def ReadCombinedPos(Small=True):
	if Small:
		fname = Globals.OutputPath + 'Venus/VenusPosSmall.bin'
	else:
		fname = Globals.OutputPath + 'Venus/VenusPos.bin'
	return RT.ReadRecarray(fname,dtype)
	
def CombinePosSmall():
	fname = Globals.OutputPath + 'Venus/VenusPosSmall.bin'
	data = ReadCombinedPos(False)
	ud = np.unique(data.Date)
	ind = np.arange(ud.size) * 24
	RT.SaveRecarray(data[ind],fname)
