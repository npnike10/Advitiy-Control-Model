import numpy as np 

def rk4(y,f,h,t):
	k1=h*f(t,y)
	k2=h*f(t+h/2,y+k1/2)
	k3=h*f(t+h/2,y+k2/2)
	k4=h*f(t+h,y+k3)

	ynext=y+1/6*(k1+2*(k2+k3)+k4)

	return ynext