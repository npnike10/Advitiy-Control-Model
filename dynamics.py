import numpy as np                                     
import numpy.linalg as la
from constants import *
#-----------------------------------------------------------------------------------------------------------------------------
def x_dot(x,control_Torque,Dist_torque):    #need m_INERTIA 
    
    #x=x.reshape((7,1))
    #control_Torque=control_Torque.reshape((3,1))
    #Dist_torque=Dist_torque.reshape((3,1))
    Tot_torque=control_Torque + Dist_torque
    q=x[0:4]; w=x[4:7]
    I_omega =np.dot(m_INERTIA,w); #I_omega=I_omega.reshape((1,3));w=w.reshape((1,3)); #m_INERTIA from constants
    m_INERTIA_inv=la.inv(m_INERTIA)
    #Tot_torque=Tot_torque.reshape((1,3))
    omega_dot=np.dot(m_INERTIA_inv,(Tot_torque - np.cross(w,I_omega)))
    #w=w.reshape((3,1))
    W=np.array([[0,w[2],-w[1],w[0]],[-w[2],0,w[0],w[1]],[w[1],-w[0],0,w[2]],[-w[0],-w[1],w[2],0]])  #notice 'W' and 'w' here 
    #q=q.reshape((4,1))
    q_dot=0.5*np.dot(W,q)    #quaternion kinematics differential equation 
    
    x = np.zeros(7)
    x[0:4] = q_dot.copy()
    
    x[4:7] = omega_dot.copy()
    
    return x   #dimension 7x1
#RK4 solver for x_dot function specially
#input x0=initial x;  h=time step;  Torques
def rk4_x(x_dot,x0,h,control_Torque,Dist_torque):
    k1 = h*x_dot(x0, control_Torque, Dist_torque)
    k2 = h*x_dot(x0+k1/2, control_Torque, Dist_torque)
    k3 = h*x_dot(x0+k2/2, control_Torque, Dist_torque)
    k4 = h*x_dot(x0+k3, control_Torque, Dist_torque)
    x1 = x0.copy()
    x1 = x1 + (k1 + 2*k2 + 2*k3 + k4)/6

    #print la.norm(x0[0:4]),la.norm(x1[0:4])
    x1[0:4] = x1[0:4]/la.norm(x1[0:4])
    
    return x1
#--------------------------------------------------------------------------------------------------------------------------------

def calculate_disturbance(sat):
    return np.array([1e-10,1e-10,1e-10])

def control_torque(control_i,v_mag_field_b):
        mag_moment_x= control_i[0]*Ax
        mag_moment_y = control_i[1]*Ay
        mag_moment_z = control_i[2]*Az
        mag_moment = mag_moment_x+mag_moment_y+mag_moment_z
       
        Torque=np.cross(mag_moment,v_mag_field_b);
        return Torque