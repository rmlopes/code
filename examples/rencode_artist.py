import ConfigParser, os
import logging
import random
import numpy as np
from numpy import array as nparray, zeros
from functools import partial
from bitstring import *
from code.utils import *
from code.utils.bitstrutils import *
from code.arnmiguel import *
import code.extendedarn as extendedarn
from code.rencode2 import *
from code.rencode import evaluatecircuit, nnlikefun, defaultnodemap
from code.utils.gui import App
from code.utils.mathlogic import *
from code.evodevo import *
from math import isnan, isinf
from code.operators import *

log = logging.getLogger(__name__)

class ArtistProb(Problem):
    funs = ['add_', 'sub_','mul_','div_',
            'sin_','cos_','sinh_', 'cosh_']
    #funs = ['gaussian', 'symmetric', 'identity','sin_']
            #'cos_','add_','cosh_','sinh_']
    terms = ['inputs[0]','inputs[1]','inputs[2]']
    arity = {}
    feedback = False
    labels = None
    def __init__(self):
        self.nout = 3
        self.ninp = 3
        Problem.__init__(self, None, defaultnodemap,
                         printmultiplecircuit)
        self.arity.update(zip(self.funs,[0]*len(self.funs)))

class TheArtist(ARNGPAgent):
        def __init__(self, config, problem, gcode = None, parentfit = 1e9):
                self.generate = extendedarn.generatechromo
                self.numproducts = 3
                ARNGPAgent.__init__(self, config, problem, gcode, parentfit)

        #def pickled(self):
         #       return self.genotype.code.bin


def render_images(pop, img_size, feedback = False, **kwargs):
        log.info('Rendering population...')
        evalf = partial(evaluatecircuit, nout = 3)
        images = []
        for i in pop:
            #log.debug('######')
            #log.debug(printcircuit(i.phenotype))
            log.debug('Rendering individual')
            striped = nparray(zeros(img_size+(3,)), dtype='int32')
            #resultdict = dict(zip([c[0] for c in i.phenotype.getcircuit()],
             #                     [1.0]*len(i.phenotype.getcircuit())))
            resultdict=dict()
            for x in range(img_size[0]):
                for y in range(img_size[1]):
                    if not feedback:
                        resultdict = dict()
                    #try:
                    r,g,b = evalf(i.phenotype.getcircuit(),
                                          nnlikefun,dict(),
                                          *prepare_inputs((x,y),img_size))

                    #print i.phenotype.getcircuit()
                    #exit(0)
                    if not isinstance(r, long) and (isnan(r) or isinf(r)):
                        r = .0
                    if abs(r) > 1.0:
                        r = r % 255
                    else:
                        r *= 255

                    if not isinstance(g, long) and (isnan(g) or isinf(g)):
                        g = .0
                    if abs(g) > 1.0:
                        g = g % 255
                    else:
                        g *= 255

                    if not isinstance(b, long) and (isnan(b) or isinf(b)):
                        b = .0
                    if abs(b) > 1.0:
                        b = b % 255
                    else:
                        b *= 255

                    try:
                            striped[x][y] = (int(abs(r)),
                                             int(abs(g)),
                                             int(abs(b)))
                    except OverflowError:
                            print "Overflow error.... "
                            striped[x][y] = 0
                    #try:
                    #if abs(r) < 1:
                     #       striped[x][y] = int(abs(r) * 255)
                    #else:
                     #       striped[x][y] = int(abs(r) % 255)
                    #except ValueError:
                     #   import traceback
                     #   print traceback.format_exc()
                     #   print 'r: ', r
                     #   exit(0)

            images.append(striped)
        log.info('done...')
        return images

def prepare_inputs(point, imgsize):
    #return normalizexy(point, imgsize) + (getcenterdistance(point,
    #imgsize),)
    center = (int(imgsize[0]/2.0),
              int(imgsize[1]/2.0))
    return point + (sqrt(pow(point[0]-center[0],2)
                         + pow(point[1]-center[1],2)),)

def getcenterdistance(point, imgsize):
    center = (int(imgsize[0]/2.0),
              int(imgsize[1]/2.0))
    return sqrt(pow(point[0]-center[0],2) + pow(point[1]-center[1],2))\
           / sqrt(pow(center[0],2) + pow(center[1],2))

def normalizexy(point, imgsize):
    return point[0]/float(imgsize[0]), point[1]/float(imgsize[1])

if __name__ == '__main__':
        #arnconfigfile = 'configfiles/arnsim.cfg'
        #log.setLevel(logging.DEBUG)
        #cfg = ConfigParser.ConfigParser()
        #cfg.readfp(open(arnconfigfile))

        p = ArtistProb()
        print p.funs
        edw = EvoDevoWorkbench(sys.argv[1],p,TheArtist)
        edw.run()

        filename = "selectedcircuits.dot"
        f = open(filename,'w')
        for i in range(5):
            f.write(edw.population[i].genotype.code.bin)
            f.write("\n")
        f.close