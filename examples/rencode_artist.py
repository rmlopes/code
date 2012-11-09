import ConfigParser, os
import logging
import random
import numpy as np
from numpy import array as nparray, zeros
from functools import partial
from bitstring import *
from code.utils import *
from code.utils.bitstrutils import *
from code.arn import *
from code.rencode import *
from code.utils.gui import App
from code.utils.mathlogic import *
from code.evodevo import *
from math import isnan, isinf
from code.operators import npoint_xover

log = logging.getLogger(__name__)

class ArtistProb(ReNCoDeProb):
    #funs = ['tanh_',
     #            'sin_','cos_','tan_']
    funs = ['add_', 'sub_', 'cos_','sin_']
    terms = ['inputs[0]','inputs[1]']
    feedback = True
    labels = None
    def __init__(self):
        ReNCoDeProb.__init__(self,None)
        #self.labels['inputs[0]'] = 'x'
        #self.labels['inputs[1]'] = 'y'
        self.bias = random.random()
#        self.funs.extend(self.extrafuns)
        self.arity.update(zip(self.funs,[0]*len(self.funs)))

def normalizexy(point, imgsize):
    return point[0]/float(imgsize[0]), point[1]/float(imgsize[1])

if __name__ == '__main__':
        #arnconfigfile = 'configfiles/arnsim.cfg'
        #log.setLevel(logging.DEBUG)
        #cfg = ConfigParser.ConfigParser()
        #cfg.readfp(open(arnconfigfile))

        p = ArtistProb()
        print p.funs
        edw = EvoDevoWorkbench(sys.argv[1],p,ReNCoDeAgent)
        edw.run()

        filename = "selectedcircuits.dot"
        f = open(filename,'w')
        #for i in edw.gui.selected:
        #f.write(printdotcircuit(edw.population[i].phenotype))
        #f.close


        #print 'Generalization: '
        #print genresult

        #proteins=[]
        #nump = 0
        #try:
        ''' f = open(sys.argv[1], 'r')
            genome = BitStream(bin=f.readline())
            arnet = ARNetwork(genome, cfg)
        except:
            #pass
            while nump < 4 or nump > 32:
                genome = BitStream(float=random.random(), length=32)
                for i in range(cfg.getint('default','initdm')):
                    genome = dm_event(genome,
                                      .02)

                    arnet = ARNetwork(genome, cfg)
                    nump = len(arnet.promlist)
        offspring = None
        popsize = 64
        pop = [ReNCoDeAgent(cfg,p) for i in range(popsize)]
        #eval_ = bindparams(cfg, evaluatecircuit)
        #plot_ = bindparams(cfg, plotindividual)

        #pop = zip()(themother, eval_(themother)),
         #      (None,0)]
        theApp = gui.App(popsize,printfun = printdotcircuit)
        theApp.images = render_images(pop, theApp, bias = p.bias)
        theApp.pop = pop
        theApp.on_execute()
        selected = []
        while theApp._running:
            if theApp.selected:
                selected = []
                selected.extend([pop[i]
                                 for i in
                                 theApp.selected[-int(sqrt(popsize)):]])
            pop = []
            for i in selected:
                pop.append(i)
                for j in range(int(sqrt(popsize))-1):
                    r = random.random()
                    offspring = ReNCoDeAgent(cfg,p,
                                             bitflipmutation(
                                                 i.genotype.code,
                                                 .01))
                    while len(offspring.genotype.promlist) < 3:
                        offspring = ReNCoDeAgent(
                            cfg, p, bitflipmutation( i.genotype.code,.01))
                    pop.append(offspring)

            theApp.images = render_images(pop, theApp, bias = p.bias)
            theApp.pop = pop
            theApp.unpause()

        filename = "selectedcircuits.dot"
        f = open(filename,'w')
        for i in selected:
            f.write(printdotcircuit(i.phenotype))
        f.close'''