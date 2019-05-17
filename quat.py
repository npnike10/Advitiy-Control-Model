import numpy as np

def cpm(v): #cross product matrix of a vector
    a1=v[0]
    a2=v[1]
    a3=v[2]
    cpm=np.array([[0,-a3,a2],
                  [a3,0,-a1],
                  [-a2,a1,0]])
    return cpm

def qinv(q): #inverse of a quaternion
    qiv=-q[1:]
    qio=q[0]

    qi=np.hstack((qio,qiv))
    return qi


def qp(q1, q2): #quaternion product
    s0=q1[0]
    s1=q1[1]
    s2=q1[2]
    s3=q1[3]

    qcross=np.array([[s0,s1,s2,s3], 
                     [s1,s0,-s3,s2],
                     [s2,s3,s0,-s1],
                     [s3,-s2,s1,s0] ])
    q=np.dot(qcross,q2) 

    return q

def rm(q): #rotation matrix corresponding to a quaternion
    q0=q[0]
    qv=q[1:]
    I=np.identity(3)
    rm=(q0**2-np.dot(qv,qv))*I+2*np.matmul(qv,qv.T)-2*q0*cpm(qv)

    return rm