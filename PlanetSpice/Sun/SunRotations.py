import numpy as np
import DateTimeTools as TT
from .Transform import HCItoIAU_SUN
from .. import Globals
from ..Tools.ContUT import ContUT
import PyFileIO as pf

dtype=[	('Date','int32'),
		('ut','float32'),
		('utc','float64')]


def SaveSolarRotations(Date0=19500101,Date1=20500101):
	'''
	List all of the times when the Sun has completed a full rotation
	
	'''
	#list the dates
	Dates = ListDates(Date0,Date1)
	nd = Dates.size	
	
	#some temporary arrays
	ut = np.zeros(nd,dtype='float32')
	x = np.zeros(nd,dtype='float32')+1.496e8
	y = np.zeros(nd,dtype='float32')
	z = np.zeros(nd,dtype='float32')
	
	#convert coordinates
	nx,ny,nz = HCItoIAU_SUN(Dates,ut,x,y,z)
	
	#calculate longitude
	lon = np.arctan2(ny,nx)*180.0/np.pi
	
	#find the crossings and interpolate
	use0 = np.where((lon[:-1] >= 0.0) & (lon[1:] < 0.0))[0]
	use1 = use0 + 1
	n = np.size(use0)
	out = np.recarray(n,dtype=dtype)
	out.Date = Dates[use0]
	m = (lon[use1]-lon[use0])/24.0
	c = lon[use0]
	out.ut = -c/m
	out.utc = ContUT(out.Date,out,ut)
	
	#save
	outpath = Globals.OutputPath + 'Sun/'
	if not os.path.isdir(outpath):
		os.system('mkdir -pv '+outpath)
	fname = outpath + '/SunRotations.dat'
	
	pf.WriteASCIIData(fname,data)

def ReadSolarRotations():
	'''
	Read the file containing the solar rotations.
	
	'''
	
	path = Globals.OutputPath + 'Sun/'
	fname = path + '/SunRotations.dat'
	
	return pf.ReadASCIIData(fname,dtype=dtype)
