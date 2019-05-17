import numpy as np

h=0.001 #constants
T=180
N=int(T/h)
kp=1
kv=1
j=np.array([[20,1.2,0.9], [1.2,17,1.4], [0.9,1.4,15]])
j1=np.linalg.inv(j) 