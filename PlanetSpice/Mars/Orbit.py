import numpy as np
from ..Tools.RotTrans import RotTrans
from ..Sun.Transform import HAEtoHCI

def OrbitHAE():
	a = 227.9392e6
	e = 0.0934
	b = a*np.sqrt(1 - e**2)
	t = np.arange(361,dtype='float32')*np.pi/180.0
	rp = (1.0 - e)*a
	x0 = a*np.cos(t)-(a-rp)
	y0 = b*np.sin(t)
	i = 5.65*np.pi/180.0
	la = 49.558*np.pi/180.0
	ap = 286.502*np.pi/180.0
	
	x1,y1 = RotTrans(x0,y0,ap)
	z1 = np.zeros(361,dtype='float32')
	
	y2,z2 = RotTrans(y1,z1,i)
	x2 = x1
	
	x3,y3 = RotTrans(x2,y2,la)
	z3 = z2
	return (x3,y3,z3)

def OrbitHCI(Date=20150101):
	x,y,z = MarsOrbitHAE()
	return HAEtoHCI(Date,0.0,x,y,z)
