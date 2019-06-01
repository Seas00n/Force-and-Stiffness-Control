from rigidBody import rigidBody
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pdb import set_trace

debug_info = 0

class armStructure(object):

    def __init__(self):
        self.TMatrixs = []#[T_{01}, T_{02}, T_{03}, ....]
        self.TMatrix = np.array([ [1,0,0,0],
                                  [0,1,0,0],
                                  [0,0,1,0],
                                  [0,0,0,1]])
        self.links = []
        self.num_links = 0
        self.joints = []
        self.jacob_matrix = np.zeros([6, 0])

    def addLink(self, joint_type, theta, d, a, alpha, link_name):
        link = rigidBody(joint_type,theta,d,a,alpha,link_name)
        self.links.append(link)
        self.num_links += 1
    '''
    kinematic 
    '''
    def resetTMatrix(self):
        self.TMatrixs = []#[T_{00},[T_{01}, T_{02}, T_{03}, ....]
        self.TMatrix = np.array([ [1,0,0,0],
                                  [0,1,0,0],
                                  [0,0,1,0],
                                  [0,0,0,1]])

    def forwardKinematics(self):
        for i in range(self.num_links):
            self.TMatrix = self.TMatrix.dot(self.links[i].AMatrix)
            if debug_info:
                print('T:      \n',self.TMatrix)
                print('A:     \n', self.links[i].AMatrix)
            self.TMatrixs.append(self.TMatrix)
            #print(self.TMatrix)
            if i==0:
                self.links[i].calculateJointCoord(np.array((0,0,0)), TMatrix=self.TMatrix)
            else:
                self.links[i].calculateJointCoord(last_joint_coord=self.links[i-1].end_joint_coord, TMatrix=self.TMatrix)

    def moveIt(self, joint_list):
        assert len(joint_list) == self.num_links-1
        for i in range(1, self.num_links):
            self.links[i].moveIt(joint_list[i-1])

        self.resetTMatrix()
        self.forwardKinematics()

    '''
    force
    '''
    def jacobianCalc(self):
        self.jacob_matrix = np.zeros([6, self.num_links-1])
        for i in range(self.num_links):
            if self.links[i].joint_type == 'nojoint':
                continue
            elif self.links[i].joint_type == 'prismatic':
                print('prismatic not yet')
                exit()
            elif self.links[i].joint_type == 'revolute':
                '''
                [j_vi] = Z_0_i-1 X (d_0_n - d_0_i-1)
                [j_wi] = Z_0_i-1
                joint i -> frame i-1
                '''      
                z_0_iminus1 = self.TMatrixs[i-1][:-1, 2]
                d_0_iminus1 = self.TMatrixs[i-1][:-1,-1]
                d_0_n       = self.TMatrix[:-1,-1]
                self.jacob_matrix[:3,i-1] = np.cross(z_0_iminus1, (d_0_n-d_0_iminus1))
                self.jacob_matrix[3:,i-1] = z_0_iminus1   
            else:
                print('error') 
                exit()

    def pseudo_inverseJacobian(self):
        return np.linalg.inv(self.jacob_matrix[:2,:])
        return self.jacob_matrix.T.dot(np.linalg.inv(self.jacob_matrix.dot(self.jacob_matrix.T)))

    def torque_displacement(self):
        '''
        F=Kx
        torque=J_transpose*F +  torque_g
        '''
        K = np.array([[1000,200],[200,2000]])
        F =np.zeros((6,26))
        delta_y = np.linspace(0, 0.025, 26).reshape((1,-1))
        delta_x = np.zeros_like(delta_y).reshape((1,-1))
        delta_displace = np.concatenate((delta_x, -delta_y), 0)
        print(delta_displace.shape)
        F[:2] = np.dot(K, delta_displace)
#        torque_g = np.array([[-4*9.8*np.cos(10/180*3.14159)*0.07],[-2*9.8*np.cos(60)*0.06]])
        torque_g = np.array([[-4*9.8*np.cos(10/180*3.14159)*0.07 - 2*9.8*(np.cos(60/180*3.14159)*0.06+np.cos(10/180*3.14159)*0.14)],
                             [            -2*9.8*np.cos(60/180*3.14159)*0.06]])

        #set_trace()
        torques = self.jacob_matrix.T.dot(F) + torque_g
#        torques = self.jacob_matrix.T.dot(F) 
        
        torques[0, np.where(torques[0]>3)]  = 3 
        torques[0, np.where(torques[0]<-3)] = -3 
        torques[1, np.where(torques[1]>2)]  = 2
        torques[1, np.where(torques[1]<-2)] = -2
        return torques,delta_y

    def visualize(self):
        x=[]
        y=[]
        z=[]
        for i in range(self.num_links):
            a,b,c = self.links[i].visualize()
            x.append(a)
            y.append(b)
            z.append(c)
        for c in range(len(x)):
            i = x[c]
            j = y[c]
            k = z[c]
            ax.plot(i,j,k,'o-')
        return x,y,z

if __name__ == '__main__':

    fig = plt.figure()
    #fig1 = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    #ax1 = fig1.add_subplot(111)

    robotic_arm = armStructure()
    robotic_arm.addLink('nojoint',    0,  0,    0,  0,'link0')
    robotic_arm.addLink('revolute',  10,  0, 0.14,  0,'link1')
    robotic_arm.addLink('revolute', 110,  0, 0.12,  0,'link2')
    #robotic_arm.addLink('revolute',   0,  10, 0,  0,'link1')
    #robotic_arm.addLink('revolute',   0,  10, 0,  -90,'link2')
    #robotic_arm.addLink('revolute',   0,  10, 0,  0,'link3')

    robotic_arm.forwardKinematics()
    x,y,z = robotic_arm.visualize()
    plt.axis('off')
    plt.pause(1)
    plt.cla()
    for i in range(20, 50, 2):
        robotic_arm.moveIt([i,110])
        x,y,z = robotic_arm.visualize()
        plt.axis('off')
        plt.pause(0.4)
        plt.cla()
    for j in range(100, 50, -2):
        robotic_arm.moveIt([50,j])
        x,y,z = robotic_arm.visualize()
        plt.axis('off')
        plt.pause(0.4)
        plt.cla()

    '''
    robotic_arm.jacobianCalc()
    torques,delta_x = robotic_arm.torque_displacement()
    lns1 = ax1.plot(delta_x[0], torques[0],   '--',  color='b',  label='tor_1a')
    lns2 = ax1.plot(delta_x[0], torques[1],  'p--',  color='b',  label='tor_2a')
    plt.yticks([i for i in range(-6,4,2)])
    plt.legend(loc='best')
    plt.ylim(-6,2.05)
    plt.grid()
    '''  
    plt.show()
