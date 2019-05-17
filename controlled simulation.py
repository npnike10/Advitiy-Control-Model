import numpy as np 
import matplotlib.pyplot as plt

def qp(q1, q2):
    b1=q1[0]
    c1=q1[1]
    d1=q1[2]
    a1=q1[3]
    b2=q2[0]
    c2=q2[1]
    d2=q2[2]
    a2=q2[3]
    q=np.array([a1*b2+b1*a2+c1*d2-d1*c2,
                a1*c2-b1*d2+c1*a2+d1*b2,
                a1*d2+b1*c2-c1*b2+d1*a2,
                a1*a2-b1*b2-c1*c2-d1*d2])
    return q

def rm(q):
    s=np.linalg.norm(q)
    q1=q[0]
    q2=q[1]
    q3=q[2]
    q4=q[3]
    rm=np.array([ [1-2*s*(q2**2+q3**2), 2*s*(q1*q2-q3*q4), 2*s*(q1*q3+q2*q4)],
        [2*s*(q1*q2+q3*q4), 1-2*s*(q1**2+q3**2), 2*s*(q2*q3-q1*q4)],
        [2*s*(q1*q3-q2*q4), 2*s*(q2*q3+q1*q4), 1-2*s*(q1**2+q2**2)]
        ])
    return rm

def fw(t):
    w=np.array([0.3*np.cos(t)*(1-np.exp(-0.01*t**2))+(0.08*np.pi+0.006*np.sin(t))*t*np.exp(-0.01*t**2),
        0.3*np.cos(t)*(1-np.exp(-0.01*t**2))+(0.08*np.pi+0.006*np.sin(t))*t*np.exp(-0.01*t**2),
        1])
    return w

def ddt(t):
    ddt=-0.3*np.sin(t)*(1-np.exp(-0.01*t**2))+0.012*np.cos(t)*t*np.exp(-0.01*t**2)+(0.08*np.pi+0.006*np.sin(t))*np.exp(-0.01*t**2)-0.02*(t**2)*np.exp(-0.01*t**2)*(0.08*np.pi+0.006*np.sin(t))
    return ddt

def wd(T):
    wd=np.array([ddt(T), ddt(T), 0])
    return wd

def cpm(v):
    a1=v[0]
    a2=v[1]
    a3=v[2]
    cpm=np.array([[0,-a3,a2],
                  [a3,0,-a1],
                  [-a2,a1,0]])
    return cpm


h=0.01 #constants
T=60
kp=1
kv=1
j=np.array([[20,1.2,0.9], [1.2,17,1.4], [0.9,1.4,15]])
j1=np.linalg.inv(j) 

t1=np.zeros(int(T/h+1))
for i in range(int(T/h)):
    t1[i+1]=t1[i]+h


qv0=np.array([0.1826,0.1826,0.1826])
qo0=-np.sqrt(1-np.linalg.norm(qv0)**2)
q=np.zeros((int(T/h+1),4))
q[0][0]=qv0[0]
q[0][1]=qv0[1]
q[0][2]=qv0[2]
q[0][3]=qo0

w0=np.array([0.1,0.1,0.1]) 
w=np.zeros((int(T/h+1),3))
w[0][0]=w0[0]
w[0][1]=w0[1]
w[0][2]=w0[2]

qd=np.zeros((int(T/h+1),4))
qd[0][0]=0
qd[0][1]=0
qd[0][2]=0
qd[0][3]=1



tau=np.zeros((int(T/h+1),3))
delw=np.zeros((int(T/h+1),3))

for i in range(int(T/h)):
    t=i*h
    qdv=qd[i][0:3]
    qdo=qd[i][3]

    qd[i+1][0:3]=qdv+h*0.5*(qdo*fw(t)+np.cross(qdv, fw(t)))
    qd[i+1][3]=qdo-h*0.5*np.dot(qdv, fw(t))

for i in range(int(T/h)):
    t=i*h
    s=qp(qd[i], q[i])
    delw1=w[i]-np.dot(rm(s), fw(t))
    delw[i]=delw1
    qv=q[i][0:3]
    qo=q[i][3]
    w1=w[i]

    q[i+1][0:3]=qv+h*0.5*(qo*w1+np.cross(qv, w1))
    q[i+1][3]=qo-h*0.5*np.dot(qv, w1)
    
    w[i+1]=w[i]+h*(-np.dot(cpm(delw1), np.dot(rm(s), fw(t)))+np.dot(rm(s), wd(t))-kp*np.dot(j1, s[0:3])-kv*np.dot(j1, delw1))

    term1=np.dot(cpm(delw1), np.dot(rm(s), fw(t)))
    tau[i]=np.cross(w[i], np.dot(j, w[i]))+np.dot(j, -term1+np.dot(rm(s), wd(t)))-kp*s[0:3]-kv*delw1

plt.figure(0)
plt.plot(t1,tau)
plt.figure(1)
plt.plot(t1,delw)