import numpy as np

def fw(t): #value of w_desired at time t
    w=np.array([0.3*np.cos(t)*(1-np.exp(-0.01*t**2))+(0.08*np.pi+0.006*np.sin(t))*t*np.exp(-0.01*t**2),
        0.3*np.cos(t)*(1-np.exp(-0.01*t**2))+(0.08*np.pi+0.006*np.sin(t))*t*np.exp(-0.01*t**2),
        1])
    return w

def wddot(t): #value of derivative of w_desired at time t
    ddt=-0.3*np.sin(t)*(1-np.exp(-0.01*t**2))+0.012*np.cos(t)*t*np.exp(-0.01*t**2)+(0.08*np.pi+0.006*np.sin(t))*np.exp(-0.01*t**2)-0.02*(t**2)*np.exp(-0.01*t**2)*(0.08*np.pi+0.006*np.sin(t))
    wddot=np.array([ddt,ddt,0])

    return wddot