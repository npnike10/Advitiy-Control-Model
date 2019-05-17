import numpy as np
import matplotlib.pyplot as plt 

#angular velocity with time

h=0.01 #time step
T=60 #total time of simulation
j=np.array([[20,1.2,0.9], [1.2,17,1.4], [0.9,1.4,15]]) #inertia matrix
w0=np.array([0.1,0.1,0.1]) #initial angular velocity
w=np.zeros((int(T/h+1),3))
w[0][0]=w0[0]
w[0][1]=w0[1]
w[0][2]=w0[2]
j1=np.linalg.inv(j)  

for i in range(int(T/h)): #propagating angular velocity using Euler's method
    wt=np.zeros((3))
    wt[0]=w[i][0]
    wt[1]=w[i][1]
    wt[2]=w[i][2]
    jw=np.matmul(j, wt.T)
    jwxw=np.cross(jw, wt)
    w[i+1][0]=(wt+h*np.matmul(j1, jwxw.T))[0]
    w[i+1][1]=(wt+h*np.matmul(j1, jwxw.T))[1]
    w[i+1][2]=(wt+h*np.matmul(j1, jwxw.T))[2]

t=np.zeros(int(T/h+1)) #defining time array
for i in range(int(T/h)):
    t[i+1]=t[i]+h

wn=np.zeros(int(T/h+1)) #defining array of norms of angular velocities at different instants
for i in range(int(T/h+1)):
    wn[i]=np.linalg.norm(w[i])

plt.figure(0)
plt.plot(t, wn)

#norm of vector part of quaternion with time

q=np.zeros((int(T/h+1), 4))
qv0=[0.1826, 0.1826, 0.1826]
qo0=-np.sqrt(1-np.linalg.norm(qv0)**2)
q[0][0]=qv0[0]
q[0][1]=qv0[1]
q[0][2]=qv0[2]
q[0][3]=qo0
w1=np.zeros(3)
qvt=np.zeros(3)

for i in range(int(T/h)):
    qvt[0]=q[i][0]
    qvt[1]=q[i][1]
    qvt[2]=q[i][2]
    qot=q[i][3]
    w1[0]=w[i][0]
    w1[1]=w[i][1]
    w1[2]=w[i][2]
    
    
    qv=qvt+h*0.5*(qot*w1+np.cross(qvt, w1))
    qo=qot-h*0.5*np.dot(qvt, w1)

    q[i+1][0]=qv[0]
    q[i+1][1]=qv[1]
    q[i+1][2]=qv[2]
    q[i+1][3]=qo
    q[i+1]=q[i+1]/np.sqrt(np.linalg.norm(qv)**2+qo**2)

qo=np.zeros(int(T/h+1))
for i in range(int(T/h+1)):
    qo[i]=q[i][3]

qvn=np.zeros(int(T/h+1))
for i in range(int(T/h+1)):
    qvn[i]=np.sqrt(q[i][0]**2+q[i][1]**2+q[i][2]**2)
    
qn=np.zeros(int(T/h+1))
for i in range(int(T/h+1)):
    qn[i]=np.linalg.norm(q[i])

plt.figure(1)
plt.plot(t, qvn)

plt.figure(2)
plt.plot(t, qo)

plt.figure(3)
plt.plot(t,qn)

plt.figure(4)
plt.plot(t, w)