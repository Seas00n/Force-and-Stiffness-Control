import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pdb import set_trace
JOINT_TYPE = ['revolute', 'prismatic', 'nojoint']
show_axis = 1

class rigidBody(object):
    '''
    link length a
    link offset d
    link twist  alpha
    joint angle theta
    '''
    #fig = plt.figure()
    #ax = Axes3D(fig)

    def __init__(self, joint_type, theta, d, a, alpha, link_name, fix_frame_coord=np.zeros([3,1])):
        self.theta = theta/180.0*np.pi
        self.d     = d
        self.a     = a
        self.alpha = alpha/180*np.pi
        self.AMatrix = self.generateHomogeneousMatrix()
        self.head_joint_coord = fix_frame_coord#x0,y0,z0
        self.end_joint_coord = np.zeros([3,1]) #x1,y1,z1
        self.joint_type = joint_type
        self.link_name = link_name
        #rigidBody.ax.set_xlabel('x')
        #rigidBody.ax.set_ylabel('y')
        #rigidBody.ax.set_zlabel('z')

    def setParameter(self, arg):
        if self.joint_type ==JOINT_TYPE[0]:
            #revolute
            self.theta = arg
        elif self.joint_type == JOINT_TYPE[1]:
            #prismatic
            self.d = arg
        else:
            #base link0
            pass

    def moveIt(self):
        pass

    def generateHomogeneousMatrix(self):
        '''
        rot_theta, trans_d: along axis Z{i-1}
        trans_a, rot_alpha: along axis X{i}
        ''' 
        rot_theta = np.array([[np.cos(self.theta), -np.sin(self.theta), 0, 0],
                              [np.sin(self.theta),  np.cos(self.theta), 0, 0],
                              [0                 ,                   0, 1, 0],
                              [0                 ,                   0, 0, 1]])

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

        #AMatrix = ((rot_theta.dot(trans_d)).dot(trans_a)).dot(rot_alpha)
        AMatrix = np.dot(rot_theta, np.dot(trans_d, np.dot(trans_a, rot_alpha)))
        AMatrix[np.where((AMatrix<0.00001) & (AMatrix>0))]=0
        AMatrix[np.where((AMatrix>-0.00001) & (AMatrix<0))]=0

        return  AMatrix

    def calculateJointCoord(self, last_joint_coord, TMatrix):
        #self.joint_coord = np.dot(self.AMatrix, last_joint_coord)
        self.head_joint_coord = last_joint_coord.reshape((3,1))
        self.end_joint_coord  = TMatrix[:3,-1].reshape((3,1))
        #set_trace()

    def visualize(self):
        #self.calculateJointCoord(last_joint_coord)
        #print(self.head_joint_coord[0], self.end_joint_coord[0])
        x = np.linspace(self.head_joint_coord[0], self.end_joint_coord[0], 100)
        y = np.linspace(self.head_joint_coord[1], self.end_joint_coord[1], 100)
        z = np.linspace(self.head_joint_coord[2], self.end_joint_coord[2], 100)
        return x,y,z
        #rigidBody.ax.plot(x,y,z,'o-')
    '''
    def show_axis(self):
        if show_axis:
            x_axis = self.AMatrix[:,0]
            y_axis = self.AMatrix[:,1]
            z_axis = self.AMatrix[:,2]
            for axis in [x_axis, y_axis, z_axis]:
                x = np.linspace(-50, 50, 100)
                y = np.linspace(self.head_joint_coord[1], self.end_joint_coord[1], 100)
                z = np.linspace(self.head_joint_coord[2], self.end_joint_coord[2], 100)
    '''
    def getInfo(self):
        print(self.link_name,'{')
        print('  theta:', self.theta, '\n  d:', self.d, '\n  a:', self.a, '\n  alpha:', self.alpha, '\n  AMatrix:', self.AMatrix)
        print('  head coord: ', self.head_joint_coord)
        print('  end coord: ', self.end_joint_coord)
        print('}')

if __name__ == '__main__':

    link0 = rigidBody('nojoint', theta=0, d=50, a=0, alpha=0, link_name='link0',head_joint_coord=np.array((1,2,3)))
    link0.visualize(np.array([0,0,0]))
    link0.getInfo()

    #link1 = rigidBody('revolute', theta=0, d=20, a=0, alpha=0, link_name='link1')
    #link1.visualize(np.array([0,0,0]))
    #link1.getInfo()
    plt.show()



