import numpy as np
from ..Tools.RotTrans import RotTrans
from ..Sun.Transform import HAEtoHCI

def OrbitHAE():
	'''
	This function returns the orbit of Venus in HAE coordinates
	(Heliocentric Aries Ecliptic) where Z is perpendicular to the 
	Earth's ecliptic plane (positive northwards), X points towards the
	first point of Aries (defined by the intersection between Earth's
	equatorial plane and the ecliptic).
	
	'''
	
	a = 108208000.0
	e = 0.006772
	b = a*np.sqrt(1 - e**2)
	t = np.arange(361,dtype='float32')*np.pi/180.0
	rp = (1.0 - e)*a
	x0 = a*np.cos(t)-(a-rp)
	y0 = b*np.sin(t)
	i = 3.86*np.pi/180.0
	la = 76.68*np.pi/180.0
	ap = 54.884*np.pi/180.0
	
	x1,y1 = RotTrans(x0,y0,ap)
	z1 = np.zeros(361,dtype='float32')
	
	y2,z2 = RotTrans(y1,z1,i)
	x2 = x1
	
	x3,y3 = RotTrans(x2,y2,la)
	z3 = z2
	return (x3,y3,z3)
	
def OrbitHCI(Date=20150101):
	'''
	Calculate the orbit positions in HCI coordinates (Heliocentric
	Inertial), where Z is the Sun's rotational axis, X is the solar 
	ascending node on the ecliptic.
	
	Inputs
	======
	Date : int
		Date in format yyyymmdd
	
	'''
	x,y,z = OrbitHAE()
	return HAEtoHCI(Date,0.0,x,y,z)
