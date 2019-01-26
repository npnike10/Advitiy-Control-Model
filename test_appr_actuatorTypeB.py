import satellite as sat
import numpy as np
import appr_actuator as act
from qnv import quatInv, quatRotate, qBI2qBO
import unittest	
from ddt import ddt,file_data,unpack

@ddt
class TestmagMoment(unittest.TestCase):
        @file_data("test-data/test_appr_actuatorTypeB.json")
        @unpack
        
        def test_Bdot(self,value):
               
                torque_req=np.asarray(value[0])
                magfield_b=np.asarray(value[1])
                app_torque_exp=np.asarray(value[2])

                vel=np.array([5,6,7])
                pos=np.array([3,4,5])

                q_BI=np.array([1,2,3])
                q_IB=quatInv(q_BI)
                q_BO=qBI2qBO(q_BI,vel,pos)

                magfield_i=quatRotate(q_IB, magfield_b)

                w_BO_b=np.array([8,6,7])
                state=np.hstack(q_BO,w_BOB)
                t0=0
                s=sat.satellite(state,t0)
                s.setMagmomentRequired_b(mag_moment_b)
                s.setMag_i(magfield_i)
                s.setControl_b(torque_req)
                
                app_torque=act.actuatorTypeB(s)

                self.assertTrue(np.allclose(app_torque, app_torque_exp))
                
               
if __name__=='__main__':
        unittest.main(verbosity=2)


