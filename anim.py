import pygame
import numpy as np
import sys
import armStructure

class Simulator(object):

	white = (255, 255, 255)

	def __init__(self):
		pygame.init()
		pygame.display.set_caption("Move")
		self.display = pygame.display.set_mode((640,480))
		self.fpsClock = pygame.time.Clock()


	def main(self):
		image = pygame.image.load("pic.png")
		radians = 0

		while 1:
			self.display.fill((255, 255, 255))

			radians += .1

			#rotated_image = pygame.transform.rotate(image, radians)
			rotated_image = pygame.transform.rotozoom(image, np.degrees(radians), 1)

			rect = rotated_image.get_rect()
			rect.center = (150, 150)

			self.display.blit(rotated_image,rect)

			for event in pygame.event.get():
				if event.type == pygame.QUIT or (event.type == pygame.KEYDOWN and event.key == pygame.K_ESCAPE):
					pygame.quit()
					sys.exit()

			pygame.display.update()
			self.fpsClock.tick(30)



if __name__ == '__main__':

	sim = Simulator()
	sim.main()

