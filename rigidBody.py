import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

JOINT_TYPE = ['revolute', 'prismatic']

class rigidBody(object):
    '''
    link length a
    link offset d
    link twist  alpha
    joint angle theta
    '''
    fig = plt.figure()
    ax = Axes3D(fig)

    def __init__(self, joint_type, theta, d, a, alpha, link_name, if_base=0):
        self.theta = theta
        self.d     = d
        self.a     = a
        self.alpha = alpha
        self.AMatrix = self.generateHomogeneousMatrix()
        self.head_joint_coord = np.zeros([3,1]) #x0,y0,z0
        self.end_joint_coord = np.zeros([3,1]) #x0,y0,z0
        self.joint_type = joint_type
        self.if_base = if_base
        self.link_name = link_name
        self.generateHomogeneousMatrix()

    def setParameter(self, arg):
        if self.joint_type ==JOINT_TYPE[0]:
            #revolute
            self.theta = arg
        else:
            #prismatic
            self.d = arg

    def moveIt(self):
        pass

    def generateHomogeneousMatrix(self):
        '''
        rot_theta, trans_d: along axis Z{i-1}
        trans_a, rot_alpha: along axis X{i}
        '''
        rot_theta = np.array([[np.cos(self.theta), -np.sin(self.theta), 0, 0],
                              [np.sin(self.theta), np.cos(self.theta),  0, 0],
                              [0,0,1,0],
                              [0,0,0,1]])

        trans_d = np.array([[1, 0, 0, 0],
                            [0, 1, 0, 0],
                            [0, 0, 1, self.d],
                            [0, 0, 0, 1]])

        trans_a = np.array([[1, 0, 0, self.a],
                            [0, 1, 0, 0],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1]])

        rot_alpha = np.array([[1, 0,                  0,                   0],
                              [0, np.cos(self.alpha), -np.sin(self.alpha), 0],
                              [0, np.sin(self.alpha), np.cos(self.alpha),  0],
                              [0, 0,                  0,                   1]])

        AMatrix = rot_theta.dot(trans_d).dot(trans_a).dot(rot_alpha)
        return  AMatrix

    def calculateJointCoord(self, last_joint_coord):
        #self.joint_coord = np.dot(self.AMatrix, last_joint_coord)
        #self.head_joint_coord = self.AMatrix[:3,-1] + last_joint_coord
        #self.end_joint_coord  = self.AMatrix[:3,:3].dot(self.head_joint_coord)*self.a
        self.head_joint_coord = last_joint_coord
        self.end_joint_coord  = self.AMatrix[:3,-1] + last_joint_coord

    def visualize(self, last_joint_coord):
        self.calculateJointCoord(last_joint_coord)
        rigidBody.ax.set_xlabel('x')
        rigidBody.ax.set_ylabel('y')
        rigidBody.ax.set_zlabel('z')
        x = np.linspace(self.head_joint_coord[0], self.end_joint_coord[0], 100)
        y = np.linspace(self.head_joint_coord[1], self.end_joint_coord[1], 100)
        z = np.linspace(self.head_joint_coord[2], self.end_joint_coord[2], 100)
        rigidBody.ax.plot(x,y,z,'o-')
        #plt.show()

    def getInfo(self):
        print(self.link_name,'{')
        print('  theta:', self.theta, '\n  d:', self.d, '\n  a:', self.a, '\n  alpha:', self.alpha, '\n  AMatrix:', self.AMatrix)
        print('  head coord: ', self.head_joint_coord)
        print('  end coord: ', self.end_joint_coord)
        print('}')

if __name__ == '__main__':
    link0 = rigidBody('revolute', theta=0, d=0, a=0, alpha=0, if_base=1, link_name='link0')
    link0.visualize(np.array([0,0,0]))
    link0.getInfo()

    link1 = rigidBody('revolute', theta=0, d=0, a=100, alpha=0, link_name='link1')
    link1.visualize(link0.end_joint_coord)
    link1.getInfo()
    plt.show()