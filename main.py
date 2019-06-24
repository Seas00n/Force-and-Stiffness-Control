from armStructure import *
from math import sqrt as sqrt
from pdb import set_trace
def build_arm():
    robotic_arm = armStructure('1')
    robotic_arm.addLink('nojoint',    0,  0,    0,  0,'link0')
    robotic_arm.addLink('revolute',  10,  0, 0.14,  0,'link1')
    robotic_arm.addLink('revolute', 110,  0, 0.12,  0,'link2')
    robotic_arm2 = armStructure('2')
    robotic_arm2.addLink('nojoint',    0,  0,    0,  0,'link0')
    robotic_arm2.addLink('revolute',  31,  0, 0.14,  0,'link1')
    robotic_arm2.addLink('revolute', 59,  0, 0.12,   0,'link2')

    robotic_arm.forwardKinematics()
    robotic_arm.jacobianCalc()
    robotic_arm2.forwardKinematics()
    robotic_arm2.jacobianCalc()
    return [robotic_arm2]

	#return [robotic_arm, robotic_arm2]

def deltdispcement_gen():
    unit = 0.001
    angle = np.linspace(0, 360/180.0*np.pi, 360)
    dd=[(unit*np.cos(theta),unit*np.sin(theta)) for theta in angle]
    #dd = [(0,-unit)]
    return dd

def  force_by_displacement(arm, dd):
        torque,F,F_xy = arm.torque_displacement(dd)
        torque = torque.flatten()
        return torque,F,F_xy


def main0():
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(211)
    ax2 = fig1.add_subplot(212)
    F_res_list=[]
    delty_list=[]
    iter_num = 26
    torques=np.zeros((2, iter_num))

    unit = 0.001  #later delete
    arms_list = build_arm()
    deltdisp_list = deltdispcement_gen()

    lns = []
    for arm in arms_list:
        for dd in deltdisp_list:
            for cnt in range(iter_num):
                torque,F, F_xy = force_by_displacement(arm, (dd[0]*cnt, dd[1]*cnt))
                update_q_list = arm.inverseKinematic(dd)
                arm.moveIt(update_q_list)
                arm.jacobianCalc()

                torques[:,cnt] = torque
                F_res_list.append(F)
                delty_list.append(-dd[1]*cnt)
            if arm.armName=='1':
                colour = 'r'
                label1 = 'tor_1a'
                label2 = 'tor_2a'
                label3 = 'A'
            else:
                colour = 'b'
                label1 = 'tor_1b'
                label2 = 'tor_2b'
                label3 = 'B'
            lns.append(ax1.plot(delty_list, torques[0],   '-',  color=colour,  label=label1))
            lns.append(ax1.plot(delty_list, torques[1],  'p-',  color=colour,  label=label2))
            lns.append(ax2.plot(delty_list, F_res_list,  'p--',  color=colour,  label=label3))

            torques=np.zeros((2, iter_num))
            F_res_list=[]
            delty_list=[]

    ax1.set_ylabel('torque')
    ax2.set_ylabel('F')
    plt.legend(loc='best')
    #plt.ylim(-6,2.05)
    plt.grid()
    plt.show()

def main1():
    fig1 = plt.figure()
    #ax1 = fig1.add_subplot(211)
    #ax2 = fig1.add_subplot(212)
    iter_num = 26
    torques=np.zeros((2, iter_num))

    arms_list = build_arm()
    deltdisp_list = deltdispcement_gen()

    lns = []

    for arm in arms_list:
        for dd in deltdisp_list:
            F_x = [];F_y=[]
            torques=np.zeros((2, iter_num))
            for cnt in range(iter_num):
                torque,F, F_xy = force_by_displacement(arm, (dd[0]*cnt, dd[1]*cnt))
                update_q_list = arm.inverseKinematic(dd)
                arm.moveIt(update_q_list)
                arm.jacobianCalc()

                torques[:,cnt] = torque
                F_x.append(F_xy[0])
                F_y.append(F_xy[1])

            plt.plot(F_x, F_y,'o')

    #plt.set_ylabel('dF_y')
    #plt.set_xlabel('dF_x')
    #plt.legend(loc='best')
    plt.xlim(-50,50)
    plt.ylim(-50,50)

    #plt.grid()
    plt.show()

if __name__ == '__main__':
    main1()

	#while True:
	#	arm.render()

'''
def deltdispcement_gen():
    unit = 0.001
    dd = [(unit,0),  (unit/sqrt(2),unit/sqrt(2)), (0,unit), (-unit/sqrt(2),unit/sqrt(2)),
          (-unit,0), (-unit/sqrt(2), -unit/sqrt(2)), (0,-unit), (unit/sqrt(2),-unit/sqrt(2)),
          (unit*sqrt(3)/2, unit/2), (-unit*sqrt(3)/2, unit/2), (-unit*sqrt(3)/2, -unit/2), (unit*sqrt(3)/2, -unit/2),
          (unit/2,unit*sqrt(3)/2), (-unit/2,unit*sqrt(3)/2),(-unit/2,-unit*sqrt(3)/2),(unit/2,-unit*sqrt(3)/2),
          (unit/2,unit*sqrt(3)/2), (-unit/2,unit*sqrt(3)/2),(-unit/2,-unit*sqrt(3)/2),(unit/2,-unit*sqrt(3)/2),
          (unit*0.93969,unit*0.342), (-unit*0.93969,unit*0.342),(-unit*0.93969,-unit*0.342),(unit*0.93969,-unit*0.342), #20
          (unit*0.9848,unit*0.342), (unit*0.9848,-unit*0.342),(-unit*0.9848,-unit*0.342),(unit*0.9848,-unit*0.342)  #10
          ]
    #dd = [(0,-unit)]
    return dd
'''
