import logging
import pygame
from pygame.locals import *
from numpy import array as nparray, zeros,ones
from pygame.surfarray import make_surface, blit_array
import random
from math import sqrt, isnan, isinf
from pygame import key
from pydot import graph_from_dot_data
from ..rencode import evaluatecircuit, nnlikefun

log = logging.getLogger(__name__)

def dumpcircuit(ind, printfun):
        circuit = ind.phenotype
        arnet = ind.genotype
        g = graph_from_dot_data(printfun(circuit,arnet = arnet))
        g.write_png('temp.png')

def zoomsurface(surface, panel):
        zoomed = pygame.transform.scale2x(surface)
        panel.image = zoomed

class App:
    def __init__(self,numimgs = 9, imgsize=(128,128),renderer = None, **kwargs):
        self.print_ = None
        try:
            self.print_ = kwargs['printfun']
        except KeyError:
                pass
        self._running = False
        self._display_surf = None
        self.size = self.weight, self.height = (int(sqrt(numimgs))*imgsize[0],
                                                int(sqrt(numimgs))*imgsize[1])
        self.gridsize = sqrt(numimgs)
        self.img_size = imgsize#(int(self.weight / self.gridsize),
                         #int(self.height / self.gridsize))
        striped = zeros((self.img_size[0], self.img_size[1],3))

        self.images = []
        #for i in range(numimgs):
         #   glevel = random.randint(0,255)
          #  striped[:] = glevel
           # self.images.append( nparray(striped, dtype='int32') )
        mainmod = __import__('__main__')
        self._render_images = getattr(mainmod, 'render_images')

        self.surfaces = []
        self.selected = []
        self.pop = []
        self.pause = False
        self.zoompanel = pygame.sprite.Sprite()
        self.zoompanel.rect = pygame.Rect(0,0,2*self.size[0],2*self.size[1])
        try:
            self.problem = kwargs['problem']
            self.print_ = self.problem.print_
        except KeyError:
            pass

    def on_init(self):
        pygame.init()
        pygame.key.set_repeat()
        pygame.display.set_caption('''Mouse:: Right - dump circuit (temp.png) :: Left - select :: space to continue''')
        self._display_surf = pygame.display.set_mode(self.size)
        self._running = True
        for i in range(len(self.pop)):
            self.surfaces.append(self._display_surf.subsurface(
                self.getXY(i)+self.img_size))

        self.images = self._render_images(self.pop,self.img_size,
                                    self.problem.feedback)
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
                dumpcircuit(self.pop[self.getIndex(event.pos)],
                            self.print_)
            else:
                print event
        if event.type == pygame.KEYUP:
            if event.key == pygame.K_SPACE:
                self.pause = True
            elif event.key == pygame.K_ESCAPE:
                self._running = False
            elif event.key == pygame.K_z:
                self.blit_images()
                pygame.display.update()
            else:
                print event
                #pause app and get new pop

        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_z:
                pos = pygame.mouse.get_pos()
                zoomsurface(self.surfaces[self.getIndex(pos)], self.zoompanel)
                self._display_surf.blit(self.zoompanel.image,
                                        self.zoompanel.rect)
                pygame.display.update()

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
            self.images = self._render_images(self.pop,self.img_size,
                                        self.problem.feedback)
            self.blit_images()

        while self._running and not self.pause:
                for event in pygame.event.get():
                    self.on_event(event)
                    self.on_loop()
                    self.on_render()

        if not self._running:
            self.on_cleanup()


class Visualizer(App):
    def __init__(self, img):
        self.img = img
        self.size = img.shape[:2]
        self._display_surf = None
        self._running = False
        self.pause = False

    def on_init(self):
        pygame.init()
        self._display_surf = pygame.display.set_mode(self.size)
        self._running = True
        blit_array(self._display_surf, self.img)
        pygame.display.flip()

if __name__ == "__main__" :
    #TODO: update this test
        theApp = App()
        theApp.on_execute()
        print theApp.selected
        theApp.unpause()
