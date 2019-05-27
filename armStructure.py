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

    def visualize(self):
        x=[]
        y=[]
        z=[]
        for i in range(self.num_links):
            a,b,c = self.links[i].visualize()
            x.append(a)
            y.append(b)
            z.append(c)
        return x,y,z

if __name__ == '__main__':

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    robotic_arm = armStructure()
    robotic_arm.addLink('nojoint',    0,  0,    0,  0,'link0')
    robotic_arm.addLink('revolute',  10,  0, 0.14,  0,'link1')
    robotic_arm.addLink('revolute', 110,  0, 0.12,  0,'link2')
    #robotic_arm.addLink('revolute',   0,  10, 0,  0,'link1')
    #robotic_arm.addLink('revolute',   0,  10, 0,  -90,'link2')
    #robotic_arm.addLink('revolute',   0,  10, 0,  0,'link3')

    robotic_arm.forwardKinematics()
    robotic_arm.jacobianCalc()
    x,y,z = robotic_arm.visualize()
    #print(robotic_arm.jacob_matrix)
    for c in range(len(x)):
        i = x[c]
        j = y[c]
        k = z[c]
        ax.plot(i,j,k,'o-')
    plt.show()

