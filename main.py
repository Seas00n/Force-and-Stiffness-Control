from armStructure import *

if __name__ == '__main__':
	arm = armStructure()
	arm.addLink('nojoint',    0,  0,    0,  0,'link0')
	arm.addLink('revolute',  10,  0, 140,  0,'link1')
	arm.addLink('revolute', 110,  0, 120,  0,'link2')
	arm.forwardKinematics()

	while True:
		arm.render()	