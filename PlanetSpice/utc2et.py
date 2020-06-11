import numpy as np
import spiceypy as sp
import DateTimeTools as TT
import os
from . import Globals

lsk_path = Globals.SpicePath + '/lsk/naif0010.tls'

def utc2et(Date,ut):
	'''
	Convert Date and ut to the ephemeris time used for SPICE.
	
	Inputs
	======
	Date : int
		Date(s) in format yyyymmdd
	ut : float
		Time(s) in hours from beginning of the day
		
	Returns
	=======
	et : ephemeris times
	
	'''
	
	#split up the dates and times
	n = np.size(ut)
	if np.size(Date) == 1:
		Date = np.zeros(n,dtype='int32') + Date
	if np.size(ut) == 1:
		ut = np.zeros(n,dtype='float32') + ut
	yr = Date//10000
	mn = (Date % 10000)//100
	dy = Date % 100
	hh,mm,ss,ms = TT.DectoHHMM(ut,ss=True,ms=True,Split=True)
	ss = np.float32(ss) + np.float32(ms)/1000
	
	#create an array of strings
	strfmt = '{:04d} {:02d} {:02d} {:02d} {:02d} {:06.3f}'
	utc_str = np.array([strfmt.format(int(yr[i]),int(mn[i]),int(dy[i]),int(hh[i]),int(mm[i]),float(ss[i])) for i in range(0,n)])

	#check that lsk is loaded
	cnt=sp.ktotal('ALL')
	loaded = False
	if cnt != 0:
		for i in range(0,cnt):
			k_name,k_type,k_src,k_handle = sp.kdata(i,'ALL',128,32,128)
			if k_name == lsk_path:
				loaded = True
				break
	
	if loaded == False:
		sp.furnsh(lsk_path)
	
	#create the output array
	et = np.zeros((n,),dtype='float64')
	for i in range(0,n):
		et[i] = sp.str2et(utc_str[i])
	
	if loaded == False:
		sp.unload(lsk_path)

	if n == 1:
		return et[0]
	else:
		return et
