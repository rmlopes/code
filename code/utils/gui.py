import pygame
from pygame.locals import *
from numpy import array as nparray, zeros,ones
from pygame.surfarray import make_surface, blit_array
import random
from math import sqrt
from pygame import key

class App:
    def __init__(self,numimgs = 9, imgsize=(128,128)):
        self._running = False
        self._display_surf = None
        self.size = self.weight, self.height = (int(sqrt(numimgs))*imgsize[0],
                                                int(sqrt(numimgs))*imgsize[1])
        self.gridsize = sqrt(numimgs)
        self.img_size = imgsize#(int(self.weight / self.gridsize),
                         #int(self.height / self.gridsize))
        striped = zeros((self.img_size[0], self.img_size[1],3))
        self.images=[]
        for i in range(numimgs):
            glevel = random.randint(0,255)
            striped[:] = glevel
            self.images.append( nparray(striped, dtype='int32') )

        self.surfaces = []
        self.selected = []
        self.pop = []
        self.pause = False

    def on_init(self):
        pygame.init()
        pygame.display.set_caption('Phenotype Selector')
        self._display_surf = pygame.display.set_mode(self.size)
        self._running = True
        for i in range(len(self.images)):
            self.surfaces.append(self._display_surf.subsurface(
                self.getXY(i)+self.img_size))

        self.blit_images()

    def blit_images(self):
        for i in range(len(self.images)):
            blit_array(self.surfaces[i], self.images[i])
        pygame.display.flip()

    def propagateselection(self):
        for i in range(len(self.images)):
            blit_array(self._display_surf.subsurface(
                self.getXY(i)+self.img_size), self.images[self.selected[-1]])
            #self._display_surf.blit(s0,(0,0))
        pygame.display.flip()

    def getXY(self, i):
        return ((i%self.gridsize)*self.img_size[0],
                int(float(i)/self.gridsize)*self.img_size[1])

    def getIndex(self, pos):
        return int(int(pos[0]/self.img_size[0]) +
                   int(pos[1]/self.img_size[1] * self.gridsize))

    def on_event(self, event):
        if event.type == pygame.MOUSEBUTTONUP:
            if event.button == 1:
                self.selected.append(self.getIndex(event.pos))
                print self.selected
            elif event.button == 3:
                print event
            else:
                print event
        if event.type == pygame.KEYUP:
            if event.key == pygame.K_SPACE:
                self.pause = True
            elif event.key == pygame.K_ESCAPE:
                self._running = False
            else:
                print event
                #pause app and get new pop
        if event.type == pygame.QUIT:
            self._running = False

    def on_loop(self):
        pass

    def on_render(self):
        pass

    def on_cleanup(self):
        pygame.quit()

    def unpause(self):
        self.pause = False
        self.selected = []
        self.on_execute()

    def on_execute(self):
        if not self._running:
            if self.on_init() == False:
                self._running = False
        else:
            self.blit_images()

        while self._running and not self.pause:
                for event in pygame.event.get():
                    self.on_event(event)
                    self.on_loop()
                    self.on_render()

        if not self._running:
            self.on_cleanup()

if __name__ == "__main__" :
        theApp = App()
        theApp.on_execute()
        print theApp.selected
        theApp.unpause()
