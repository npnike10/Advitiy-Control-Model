import numpy as np
from qdesired import qd
from constants import h,j,j1,kp,kv
from wdesired import fw,wddot
from quat import qp,qinv,rm,cpm


def ydot(t,y): #the function which is the derivative of state vector 
    qo=y[0]
    qv=y[1:4]
    w=y[4:]
    q_desired=qd(int(t/h))
    s=qp(qinv(q_desired),y[0:4])
    delw=w-np.dot(rm(s),fw(t))

    f=np.zeros(7)

    f[0]=-0.5*np.dot(qv,w)
    f[1:4]=0.5*(qo*w+np.cross(qv,w))
    f[4:]=-np.dot(cpm(delw),np.dot(rm(s),fw(t)))+np.dot(rm(s),wddot(t))-kp*np.dot(j1,s[1:])-kv*np.dot(j1,delw)

    return f

def y2dot(t,y):
    so=y[0]
    sv=y[1:4]
    dw=y[4:]

    f=np.zeros(7)

    f[0]=-0.5*np.dot(sv,dw)
    f[1:4]=0.5*(so*dw+np.cross(sv,dw))
    f[4:]=-kp*np.dot(j1,sv)-kv*np.dot(j1,dw)

    return f
