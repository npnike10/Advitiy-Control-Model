import numpy as np 
import matplotlib.pyplot as plt
from constants import h,T,kp,kv,j,N
from qdesired import qd
from quat import qp,qinv,rm,cpm
from wdesired import fw,wddot
import rk4
from dynamics import y2dot

t1=np.zeros(N+1)   #array of time instants
for i in range(int(T/h)):
    t1[i+1]=t1[i]+h

qv0=np.array([0.1826,0.1826,0.1826]) #initial q_actual components
qo0=-np.sqrt(1-np.linalg.norm(qv0)**2) 
q0=np.hstack((qo0,qv0))
w0=np.array([0.1,0.1,0.1]) #initial w_actual components

state=np.zeros((N+1,7)) #array of state vectors at different time instants
state[0][0:4]=qp(qinv(qd(0)),q0)
state[0][4:]=w0-np.dot(rm(state[0][0:4]),fw(0))

tau=np.zeros((N+1,3)) #array of torque applied at different time instants

for i in range(N): #propagating state
    time=i*h
    w=state[i][4:]+np.dot(rm(state[i][0:4]),fw(time))
    delw=state[i][4:]
    q_error=state[i][0:4]

    state[i+1]=rk4.rk4(state[i],y2dot,h,time)
    state[i+1][0:4]=state[i+1][0:4]/np.linalg.norm(state[i+1][0:4])

    tau[i+1]=np.cross(w,np.dot(j,w))+np.dot(j,-np.dot(cpm(delw),np.dot(rm(q_error),fw(time)))+np.dot(rm(q_error),wddot(time)))-kp*q_error[1:]-kv*delw

dw=np.zeros((N+1,3))
for i in range(N+1):
    dw[i]=state[i][4:]

print(dw[90000])
plt.figure(0)
plt.plot(t1,tau)
plt.savefig('torque_q_BO.png')
plt.figure(1)
plt.plot(t1,dw)
plt.savefig('w_BO_q_BO.png')


