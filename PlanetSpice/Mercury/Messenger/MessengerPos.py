import numpy as np
import spiceypy as sp
from ..utc2et import utc2et
from scipy.interpolate import InterpolatedUnivariateSpline
import os
import DateTimeTools as TT
from ...Tools.FileSearch import FileSearch
from ...Tools.ContUT import ContUT
from ... import Globals

lsk_path = Globals.SpicePath + '/lsk/naif0010.tls'
spk_kernel = Globals.SpicePath + '/bodies/de432s.bsp'
#spk_kernel = Globals.SpicePath + '/messenger/spk/de405.bsp'
spk_kernel2 = Globals.SpicePath + '/messenger/spk/msgr_20040803_20150501_od423sc_0.bsp'
pck_kernel = Globals.SpicePath + '/bodies/pck00010.tpc'

hci_kernel = Globals.SpicePath + '/frames/tk/sunframes.tk'
mso_kernel = Globals.SpicePath + '/frames/tk/MercuryMSO.tk'
sclk_kernel = Globals.SpicePath + '/messenger/sclk/messenger_2339.tsc'
ik_kernel = Globals.SpicePath + '/messenger/ik/msgr_grns_v110.ti'
fk_kernel = Globals.SpicePath + '/messenger/fk/msgr_v231.tf'
ck_path = Globals.SpicePath + '/messenger/ck/'


	#position dtype
dtype = [	('Date','int32'),
			('ut','float32'),
			('utc','float64'),
			('x','float32'),
			('y','float32'),
			('z','float32')]

	


def MET(Date,ut):
	'''
	This might return Mission Elapsed Time, but who knows if it actually
	works!
	
	'''
	
	
	ut = np.array([ut]).flatten()
	n = ut.size
	if np.size(Date) == 1:
		Date = np.array([Date]*n)
	else:
		Date = np.array([Date]).flatten()

	et = np.zeros((n,),dtype='float64')
	
	sp.furnsh(lsk_path)
	sp.furnsh(sclk_kernel)


	ud = np.unique(Date)
	for i in range(0,ud.size):
		use = np.where(Date == ud[i])[0]
		tmp = utc2et(ud[i],0.0)
		et[use] = tmp + ut[use]*3600.0

	met = np.zeros(n,dtype='float64')
	for i in range(0,n):
		met[i] = sp.sce2t(-236,et[i])
	
	return met/1e6
	
				

def OrientationMSO(Date,ut,Verbose=False):
	'''
	This should return the direction in which MESSENGER is oriented
	
	no idea if this is right
	
	'''
	
	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
	
	#find the ck kernels
	ck = 'msgr_{:04d}_v*.bc'
	
	yymm = (Date//100) % 10000
	uyymm = np.unique(yymm)
	nck = np.size(uyymm)
	ck_kernel = []
	for i in range(0,nck):
		files = FileSearch(ck_path,ck.format(uyymm[i]))
		if files.size > 0:
			ck_kernel.append(ck_path+files[-1])
	nck = np.size(ck_kernel)
	

	#load all the kernels
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(spk_kernel2)
	sp.furnsh(sclk_kernel)
	sp.furnsh(pck_kernel)
	sp.furnsh(mso_kernel)	
	sp.furnsh(ik_kernel)
	sp.furnsh(fk_kernel)
	for i in range(0,nck):
		sp.furnsh(ck_kernel[i])

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
	if Verbose:
		for i in range(0,n):
			print('\rVector {0} of {1}'.format(i+1,n),end='')
			pos,lt = sp.spkpos('MERCURY',et[i],'MSGR_SPACECRAFT','NONE','MESSENGER')
			x[i] = pos[0]
			y[i] = pos[1]
			z[i] = pos[2]
		print()
	else:
		for i in range(0,n):
			pos,lt = sp.spkpos('MERCURY',et[i],'MSGR_SPACECRAFT','NONE','MESSENGER')
			x[i] = pos[0]
			y[i] = pos[1]
			z[i] = pos[2]	
	
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(spk_kernel2)
	sp.unload(sclk_kernel)
	sp.unload(pck_kernel)
	sp.unload(mso_kernel)
	sp.unload(ik_kernel)
	sp.unload(fk_kernel)
	for i in range(0,nck):
		sp.unload(ck_kernel[i])

	return (x,y,z)

def NSOrientationMSO(Date,ut):
	'''
	Get the orientation of MESSENGER NS?
	'''

	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
	
	
	ck = 'msgr_{:04d}_v*.bc'
	
	yymm = (Date//100) % 10000
	uyymm = np.unique(yymm)
	nck = np.size(uyymm)
	ck_kernel = []
	for i in range(0,nck):
		files = FileSearch(ck_path,ck.format(uyymm[i]))
		if files.size > 0:
			ck_kernel.append(ck_path+files[-1])
	nck = np.size(ck_kernel)
	
	#load kernels
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(spk_kernel2)
	sp.furnsh(sclk_kernel)
	sp.furnsh(pck_kernel)
	sp.furnsh(mso_kernel)	
	sp.furnsh(ik_kernel)
	sp.furnsh(fk_kernel)
	for i in range(0,nck):
		sp.furnsh(ck_kernel[i])

	if np.size(Date) == 1 :
		et[0] = utc2et(Date,0.0)
		et = et[0] + ut*3600.0
	else:
		ud = np.unique(Date)
		for i in range(0,ud.size):
			use = np.where(Date == ud[i])[0]
			tmp = utc2et(ud[i],0.0)
			et[use] = tmp + ut[use]*3600.0

	#calculate positions
	m = sp.pxform('MERCURYMSO','MSGR_GRNS_NS',et)

		
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(spk_kernel2)
	sp.unload(sclk_kernel)
	sp.unload(pck_kernel)
	sp.unload(mso_kernel)
	sp.unload(ik_kernel)
	sp.unload(fk_kernel)
	for i in range(0,nck):
		sp.unload(ck_kernel[i])

	return m

def OrientationSUN(date,ut):

	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
	
	
	ck = 'msgr_{:04d}_v*.bc'
	
	yymm = (Date//100) % 10000
	uyymm = np.unique(yymm)
	nck = np.size(uyymm)
	ck_kernel = []
	for i in range(0,nck):
		files = FileSearch(ck_path,ck.format(uyymm[i]))
		if files.size > 0:
			ck_kernel.append(ck_path+files[-1])
	nck = np.size(ck_kernel)
	
	#load kernels
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(spk_kernel2)
	sp.furnsh(sclk_kernel)
	sp.furnsh(pck_kernel)
	sp.furnsh(mso_kernel)	
	sp.furnsh(ik_kernel)
	sp.furnsh(fk_kernel)
	for i in range(0,nck):
		sp.furnsh(ck_kernel[i])

	if np.size(Date) == 1 :
		et[0] = utc2et(Date,0.0)
		et = et[0] + ut*3600.0
	else:
		ud = np.unique(Date)
		for i in range(0,ud.size):
			use = np.where(Date == ud[i])[0]
			tmp = utc2et(ud[i],0.0)
			et[use] = tmp + ut[use]*3600.0

			
	#m = sp.pxform('J2000','MSGR_SPACECRAFT',et[0])
	for i in range(0,n):
		pos,lt=sp.spkpos('SUN',et[i],'MSGR_SPACECRAFT','NONE','MESSENGER')
		x[i]=pos[0]
		y[i]=pos[1]
		z[i]=pos[2]
		
		
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(spk_kernel2)
	sp.unload(sclk_kernel)
	sp.unload(pck_kernel)
	sp.unload(mso_kernel)
	sp.unload(ik_kernel)
	sp.unload(fk_kernel)
	for i in range(0,nck):
		sp.unload(ck_kernel[i])

	return (x,y,z)


def PosMSM(Date,ut):
	'''
	Messenger position in MSM coords (km)
	
	'''
	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
		
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(spk_kernel2)
	sp.furnsh(pck_kernel)
	sp.furnsh(mso_kernel)	
	
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
		pos,lt = sp.spkpos('MESSENGER',et[i],'MERCURYMSO','NONE','MERCURY')
		x[i] = pos[0]
		y[i] = pos[1]
		z[i] = pos[2]-478.0		
	
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(spk_kernel2)
	sp.unload(pck_kernel)
	sp.unload(mso_kernel)

	return (x,y,z)


def PosHCI(Date,ut):
	'''
	Messenger position in HCI coords (km)
	
	'''	

	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
		
		
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(spk_kernel2)
	sp.furnsh(pck_kernel)
	sp.furnsh(hci_kernel)	
	
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
		pos,lt = sp.spkpos('MESSENGER',et[i],'HCI','NONE','SUN')
		x[i] = pos[0]
		y[i] = pos[1]
		z[i] = pos[2]		
	
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(spk_kernel2)
	sp.unload(pck_kernel)
	sp.unload(hci_kernel)

	return (x,y,z)

def CarringtonLongitude(Date,ut):
	'''
	Get MESSENGER's Carrington longitude
	
	'''

	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	lon=np.zeros(n,dtype='float64')

		
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(spk_kernel2)
	sp.furnsh(pck_kernel)
	
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
		pos,lt = sp.spkpos('MESSENGER',et[i],'IAU_SUN','NONE','SUN')
		lon[i] = np.arctan2(pos[1],pos[0])
	
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(spk_kernel2)
	sp.unload(pck_kernel)

	return (lon)

	

def PosHAE(date,ut):
	'''
	Messenger position in HAE coords (km)
	
	'''	

	n = np.size(ut)

	if np.size(ut) == 1:
		ut = np.array([ut])
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
		
	sp.furnsh(lsk_path)
	sp.furnsh(spk_kernel)
	sp.furnsh(spk_kernel2)
	sp.furnsh(pck_kernel)
	sp.furnsh(hci_kernel)	
	
	if np.size(Date) == 1 :
		et[0] = utc2et(Date,0.0)
		et = et[0] + ut*3600.0
	else:
		ud = np.unique(Date)
		for i in range(0,ud.size):
			use = np.where(Date == ud[i])[0]
			tmp = utc2et(ud[i],0.0)
			et[use] = tmp + ut[use]*3600.0
			
	
	pos,lt = sp.spkpos('MESSENGER',et,'ECLIPDATE','NONE','SUN')
	x = pos.T[0]
	y = pos.T[1]
	z = pos.T[2]		
	sp.unload(lsk_path)
	sp.unload(spk_kernel)
	sp.unload(spk_kernel2)
	sp.unload(pck_kernel)
	sp.unload(hci_kernel)

	return (x,y,z)

def HAELon(Date,ut):
	'''
	Messenger lingitude in HCI coords
	
	'''	

	x,y,z = MessengerPosHAE(Date,ut)
	return np.arctan2(y,x)*180.0/np.pi


def SaveMinutePos():
	path = Globals.OutputPath + 'Mercury/MESSENGER/MinutePos/'
	StartDate = 20110323
	EndDate = 20150431
	
	date = StartDate
	ut = np.arange(1440.0)/60.0
	while date <= EndDate:
		fname = path + '{:08d}.bin'.format(date)
		x,y,z = MessPosMSM(date,ut)
		
		data = np.recarray(1440,dtype=dtype)
		data.Date = Date
		data.ut = ut
		data.utc = ContUT(data.Date,data.ut)
		data.x = x/2440.0
		data.y = y/2440.0
		data.z = z/2440.0
		
		RT.SaveRecarray(data,fname)
		
		date = TT.PlusDay(date)
