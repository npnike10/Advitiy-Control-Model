import numpy as np
#import os
import matplotlib.pyplot as plt

path = os.path.abspath()

state = np.genfromtxt('state.csv', delimiter=',')
t_p = np.genfromtxt('time_p.csv', delimiter=',') 
t = []
init = t_p[0]
end = t_p[1]
h = t_p[2]
N = int(((end-init))/h)
for i in range(0,N+1):
	t.append(init+i*h)

plt.figure(1)
plt.plot(t, state[ : , 0])
plt.xlabel('time')
plt.ylabel('q1')
plt.savefig('q1 vs t')

#my_file = 'q1 vs t.png'
#plt.savefig(os.path.join(path, my_file))

plt.figure(2)
plt.plot(t, state[ : , 1])
plt.xlabel('time')
plt.ylabel('q2')
plt.savefig('q2 vs t')

plt.figure(3)
plt.plot(t, state[ : , 2])
plt.xlabel('time')
plt.ylabel('q3')
plt.savefig('q3 vs t')

plt.figure(4)
plt.plot(t, state[ : , 3])
plt.xlabel('time')
plt.ylabel('q4')
plt.savefig('q4 vs t')

plt.figure(5)
plt.plot(t, state[ : , 4])
plt.xlabel('time')
plt.ylabel('w1')
plt.savefig('w1 vs t')

plt.figure(6)
plt.plot(t, state[ : , 5])
plt.xlabel('time')
plt.ylabel('w2')
plt.savefig('w2 vs t')

plt.figure(7)
plt.plot(t, state[ : , 6])
plt.xlabel('time')
plt.ylabel('w3')
plt.savefig('w3 vs t')

plt.figure(8)
plt.plot(t, state[ : , ])
plt.xlabel('time')
plt.ylabel('w3')
plt.savefig('w3 vs t')

