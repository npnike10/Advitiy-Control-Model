import numpy as np 
import matplotlib.pyplot as plt
import rk4
from quat import qp,qinv,cpm,rm
from dynamics import ydot
from qdesired import qd
from constants import T,h,j,kp,kv
from wdesired import fw,wddot



t1=np.zeros(int(T/h+1))   #array of time instants
for i in range(int(T/h)):
    t1[i+1]=t1[i]+h


qv0=np.array([0.1826,0.1826,0.1826]) #initial q_actual components
qo0=-np.sqrt(1-np.linalg.norm(qv0)**2)
w0=np.array([0.1,0.1,0.1]) #initial w_actual components
state=np.zeros((int(T/h+1),7)) #array of state vectors at different time instants
state[0][0:4]=np.hstack((qo0,qv0))
state[0][4:]=w0


tau=np.zeros((int(T/h+1),3)) #array of torque applied at different time instants
delw=np.zeros((int(T/h+1),3)) #array of delta w at different time instants
q_error=np.zeros((int(T/h+1),4)) #array of error quaternion at different time instants 


for i in range(int(T/h)): #propagating state
    time=i*h
    q_error[i]=qp(qinv(qd(i)),state[i][0:4])
    delw[i][:]=state[i][4:]-np.dot(rm(q_error[i]),fw(time))
    
    state[i+1]=rk4.rk4(state[i],ydot,h,time)
    state[i+1][0:4]=state[i+1][0:4]/np.linalg.norm(state[i+1][0:4])

    tau[i+1]=np.cross(state[i][4:],np.dot(j,state[i][4:]))+np.dot(j,-np.dot(cpm(delw[i]),np.dot(rm(q_error[i]),fw(time)))+np.dot(rm(q_error[i]),wddot(time)))-kp*q_error[i][1:]-kv*delw[i]
    
print(delw[86000])
plt.figure(0)
plt.plot(t1,tau)
#plt.savefig('torque_q_BI.png')
plt.figure(1)
plt.plot(t1,delw)
#plt.savefig('w_BOB_q_BI.png')
