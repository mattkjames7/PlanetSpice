import numpy as np
import spiceypy as sp
from ...utc2et import utc2et
from scipy.interpolate import InterpolatedUnivariateSpline
import os
from ... import Globals
import RecarrayTools as RT
import DateTimeTools as TT
from ...Tools.ContUT import ContUT

lsk_path = Globals.SpicePath + '/bepi/misc/kernels/lsk/naif0011.tls'
pck_kernel = Globals.SpicePath + '/bodies/pck00010.tpc'
sclk_kernel = Globals.SpicePath + '/bepi/misc/kernels/sclk/bc_mpo_sclk_20150124_fake.tsc'
mmo_kernel = Globals.SpicePath + '/bepi/misc/kernels/spk/BC-ESC-DF-50028_9133MMO.bsp'
mpo_kernel = Globals.SpicePath + '/bepi/misc/kernels/spk/BC-ESC-DF-50027_9217MPO.bsp'
de430_kernel = Globals.SpicePath + '/bepi/misc/kernels/spk/de430.bsp'
de432s_kernel = Globals.SpicePath + '/bepi/misc/kernels/spk/de432s.bsp'
hci_kernel = Globals.SpicePath + '/frames/tk/sunframes.tk'
mso_kernel = Globals.SpicePath + '/frames/tk/MercuryMSO.tk'

#position dtype
dtype = [	('Date','int32'),
			('ut','float32'),
			('utc','float64'),
			('x','float32'),
			('y','float32'),
			('z','float32')]


#currently the spice kernels allow times from
#2025-03-27T19:28:24.96000004 - 2027-06-12T18:06:36.05760002 (MPO)
#2025-01-02T12:13:52.03200004 - 2027-03-19T10:59:30.59520000 (MMO)

#I will use 20250328 - 20270318
def MPOPosMSM(Date,ut):
	'''
	Position of MPO in MSM coords
	
	'''
	
	#create the output arrays
	if np.size(np.shape(ut)) == 0:
		ut = np.array([ut])
	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
		
	#load kernels
	sp.furnsh(lsk_path)
	sp.furnsh(sclk_kernel)
	sp.furnsh(de430_kernel)
	sp.furnsh(mpo_kernel)
	sp.furnsh(pck_kernel)
	sp.furnsh(mso_kernel)	
	
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
			
	#get positions
	pos,lt = sp.spkpos('MPO',et,'MERCURYMSO','NONE','MERCURY')
	x = pos.T[0]
	y = pos.T[1]
	z = pos.T[2]-478.0			
	
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(sclk_kernel)
	sp.unload(de430_kernel)
	sp.unload(mpo_kernel)
	sp.unload(pck_kernel)
	sp.unload(mso_kernel)

	return (x,y,z)
	
	
def MMOPosMSM(Date,ut):
	'''
	Position of MMO in MSM coords
	
	'''
	
	#create the output arrays
	if np.size(np.shape(ut)) == 0:
		ut = np.array([ut])
	n = np.size(ut)
	et = np.zeros((n,),dtype='float64')
	x = np.zeros(n,dtype='float64')
	y = np.zeros(n,dtype='float64')
	z = np.zeros(n,dtype='float64')
		
	#load kernels
	sp.furnsh(lsk_path)
	sp.furnsh(sclk_kernel)
	sp.furnsh(de430_kernel)
	sp.furnsh(mpo_kernel)
	sp.furnsh(pck_kernel)
	sp.furnsh(mso_kernel)	
	
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
			
	#get positions
	pos,lt = sp.spkpos('MMO',et,'MERCURYMSO','NONE','MERCURY')
	x = pos.T[0]
	y = pos.T[1]
	z = pos.T[2]-478.0			
	
	#unload kernels
	sp.unload(lsk_path)
	sp.unload(sclk_kernel)
	sp.unload(de430_kernel)
	sp.unload(mpo_kernel)
	sp.unload(pck_kernel)
	sp.unload(mso_kernel)

	return (x,y,z)
	
def SavePosMin(sc,StartDate=20250328,StopDate=20270318):
	'''
	Save the positions of MMO/MPO 
	
	'''
	
	outpath = Globals.OutputPath + 'Mercury/Bepi/pos/{:s}/min/'.format(sc.upper())
	if not os.path.isdir(outpath):
		os.system('mkdir -pv '+outpath)
	
	ut = np.arange(24.0*60.0)/60.0
	n = np.int32(np.size(ut))
	if sc.upper() == 'MMO':
		POSFUN = MMOPosMSM
	elif sc.upper() == 'MPO':
		POSFUN = MPOPosMSM
	else:
		return
	Date = StartDate
	while Date <= StopDate:
		print('\r{:08}'.format(Date),end='')
		fname = outpath + '{:08d}.bin'.format(Date)
		x,y,z = POSFUN(Date,ut)
		
		
		data = np.recarray(1440,dtype=dtype)
		data.Date = Date
		data.ut = ut
		data.utc = ContUT(data.Date,data.ut)
		data.x = x/2440.0
		data.y = y/2440.0
		data.z = z/2440.0
		
		RT.SaveRecarray(data,fname)
		
		Date = TT.PlusDay(Date)
		
	print('')


def ReadPos(Date,sc,Minute=True):
	'''
	Read the positions of MMO/MPO
	'''
	if Minute:
		path = Globals.OutputPath + 'Mercury/Bepi/pos/{:s}/min/'.format(sc.upper())
	else:
		path = Globals.OutputPath + 'Mercury/Bepi/pos/{:s}/sec/'.format(sc.upper())
	fname = path + '{:08d}.bin'.format(Date)
	return RT.ReadRecarray(fname,dtype)
	
	
def SavePosSec(sc,StartDate=20250328,StopDate=20270318):
	'''
	Save the positions of MMO/MPO 
	
	
	'''
	outpath = Globals.OutputPath + 'Mercury/Bepi/pos/{:s}/sec/'.format(sc.upper())
	if not os.path.isdir(outpath):
		os.system('mkdir -pv '+outpath)
			
	ut = np.arange(24.0*3600.0)/3600.0
	n = np.int32(np.size(ut))

	Date = StartDate
	while Date <= StopDate:
		print('\r{:08}'.format(Date),end='')
		fname = outpath + '{:08d}.bin'.format(Date)
		pos = ReadPos(Date,sc,True)

		fx = InterpolatedUnivariateSpline(pos.ut,pos.x)
		fy = InterpolatedUnivariateSpline(pos.ut,pos.y)
		fz = InterpolatedUnivariateSpline(pos.ut,pos.z)
		
		x = fx(ut)
		y = fy(ut)
		z = fz(ut)
		
		data = np.recarray(86400,dtype=dtype)
		data.Date = Date
		data.ut = ut
		data.utc = ContUT(data.Date,data.ut)
		data.x = x
		data.y = y
		data.z = z
		
		RT.SaveRecarray(data,fname)
		
		Date = TT.PlusDay(Date)
		
	print('')
