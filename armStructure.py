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
        self.kmatrix = np.array([[1000,200],[200,2000]])

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

    def inverseKinematic(self, delty):
        #if (np.linalg.matrix_rank(self.jacob_matrix) == self.jacob_matrix.shape[0]) and np.linalg.matrix_rank(self.jacob_matrix)  == self.jacob_matrix.shape[1]
        inverseJacobianXYPos = np.linalg.inv(self.jacob_matrix[:2])
        delta_joint_variable = inverseJacobianXYPos.dot(np.array([[0],[delty]]))
        q_list = []
        for i in range(self.num_links-1):
            res = self.links[i+1].getJointAngle() + float(delta_joint_variable[i])
            q_list.append(res)
        return q_list


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
        #print('jacob_matrix:   \n', self.jacob_matrix)


    def pseudo_inverseJacobian(self):
        #return np.linalg.inv(self.jacob_matrix[:2,:])
        return self.jacob_matrix.T.dot(np.linalg.inv(self.jacob_matrix.dot(self.jacob_matrix.T)))


    
    def torque_displacement(self, delta_y):
        '''
        F=Kx
        torque=J_transpose*F +  torque_g
        '''
        G1 = 4*9.8
        G2 = 2*9.8
        theta1 = self.links[1].getJointAngle()/180.0*np.pi
        theta2 = self.links[2].getJointAngle()/180.0*np.pi
        length1 = 0.07
        length2 = 0.06

        F =np.zeros((6,1))
        delta_y = delta_y
        delta_x = 0
        delta_displace = np.array([[delta_x],[delta_y]])

        F[:2] = np.dot(self.kmatrix, delta_displace)

        jacobian_com1 = np.array([[-np.sin(theta1)*length1, 0],
                                  [np.cos(theta1)*length1,  0]])
        jacobian_com2 = np.array([[-np.sin(theta1)*length1*2-np.sin(theta1+theta2)*length2, -np.sin(theta1+theta2)*length2],
                                  [ np.cos(theta1)*length1*2+np.cos(theta1+theta2)*length2,  np.cos(theta1+theta2)*length2]])
        '''
        torque_g = np.array([[-G1*np.cos(theta1)*length1 - G2* ( np.cos(theta1)*length1*2 + np.cos(theta1+theta2)*length2 )],
                             [G2*np.cos(theta1+theta2)*length2]])
        torques = self.jacob_matrix.T.dot(F) + torque_g
        '''
        
        #torque_g  = jacobian_com1.T.dot(np.array([[G1],[0]])) + jacobian_com2.T.dot(np.array([[0],[G2]]))
        torque_g  = jacobian_com1.T.dot(np.array([[0],[-G1]])) + jacobian_com2.T.dot(np.array([[0],[-G2]]))
        torques = self.jacob_matrix.T.dot(F)# + torque_g
        
        #torques = self.jacob_matrix.T.dot(F) #+ np.array([[0.9],[0.9]])
        
        torques[0, np.where(torques[0]>3)]  = 3 
        torques[0, np.where(torques[0]<-3)] = -3 
        torques[1, np.where(torques[1]>2)]  = 2
        torques[1, np.where(torques[1]<-2)] = -2
        #set_trace()
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

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    #fig = plt.figure()
    #ax = Axes3D(fig)
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
    #ax.set_zlabel('z')

    robotic_arm = armStructure()
    robotic_arm.addLink('nojoint',    0,  0,    0,  0,'link0')
    robotic_arm.addLink('revolute',  10,  0, 0.14,  0,'link1')
    robotic_arm.addLink('revolute', 110,  0, 0.12,  0,'link2')

    #### force
    iter_num = 26
    torques=np.zeros((2, iter_num))
    delt_list=[]
    delty = 0
    unit = 0.025/(iter_num-1)
    robotic_arm.forwardKinematics()
    robotic_arm.jacobianCalc()
    #set_trace()
    for i in range(26):
        torque,delta_x = robotic_arm.torque_displacement(delty)
        torques[:,i] = torque.flatten()
        delt_list.append(-delty)      
        delty += -unit
        q_list = robotic_arm.inverseKinematic(-unit)
        robotic_arm.moveIt(q_list)
        robotic_arm.jacobianCalc()

    '''
    # kinetics
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
    #set_trace()
    lns1 = ax1.plot(delt_list, torques[0],   '--',  color='b',  label='tor_1a')
    lns2 = ax1.plot(delt_list, torques[1],  'p--',  color='b',  label='tor_2a')

    '''

    robotic_arm2 = armStructure()
    robotic_arm2.addLink('nojoint',    0,  0,    0,  0,'link0')
    robotic_arm2.addLink('revolute',  31,  0, 0.14,  0,'link1')
    robotic_arm2.addLink('revolute', 59,  0, 0.12,  0,'link2')
    #### force
    iter_num = 26
    torques=np.zeros((2, iter_num))
    delt_list=[]
    delty = 0
    unit = 0.025/(iter_num-1)
    robotic_arm2.forwardKinematics()
    robotic_arm2.jacobianCalc()
    #set_trace()
    for i in range(26):
        torque,delta_x = robotic_arm2.torque_displacement(delty)
        torques[:,i] = torque.flatten()
        delt_list.append(-delty)      
        delty += -unit
        q_list = robotic_arm2.inverseKinematic(-unit)
        robotic_arm2.moveIt(q_list)
        robotic_arm2.jacobianCalc()

    lns3 = ax1.plot(delt_list, torques[0],   '-',  color='r',  label='tor_1b')
    lns4 = ax1.plot(delt_list, torques[1],  'p-',  color='r',  label='tor_2b')
    '''

#############################
    #plt.yticks([i for i in range(-6,4,2)])
    plt.legend(loc='best')
    #plt.ylim(-6,2.05)
    plt.grid()
      
    plt.show()
