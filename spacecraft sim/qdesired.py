import numpy as np
from constants import T,h
import rk4
from wdesired import fw


def qddot(t,qd):
    qdv=qd[1:]
    qdo=qd[0]
    f=np.zeros(4)
    f[0]=-0.5*np.dot(qdv,fw(t))
    f[1:]=0.5*(qdo*fw(t)+np.cross(qdv,fw(t)))

    return f

qdes=np.zeros((int(T/h+1),4)) #array of q_desired at different time instants
qdes[0][:]=np.array([1,0,0,0])

for i in range(int(T/h)): #propagating q_desired
    tn=i*h
    
    qdes[i+1][:]=rk4.rk4(qdes[i],qddot,h,tn)
    qdes[i+1][:]=qdes[i+1][:]/np.linalg.norm(qdes[i+1][:])

def qd(n):
	return qdes[n]