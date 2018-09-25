import ConfigParser, os
import logging
import random
import numpy as np
from numpy import array as nparray, zeros
from functools import partial
from bitstring import *
from code.utils import *
from code.utils.bitstrutils import *
from code.extendedarn import *
from code.cellcode import *
from code.utils import gui

def render_images(pop, app, **kwargs):
        log.info('Rendering popoulation...')
        #ind.arn.nstepsim(2000)#, *inputs)
        #get outputs
        n = 3
        ok=0
        images = []
        for i in pop:
            log.debug('Rendering individual')
            striped = nparray(zeros(app.img_size+(3,)), dtype='int32')
            for x in range(app.img_size[0]):
                for y in range(app.img_size[1]):
                    i.phenotype.reset()
                    #print 'MAX = ',app.img_size
                    i.phenotype.simulate(*normalizetocc((x,y),app.img_size))
                    striped[x][y] = getoutput(i.phenotype)
            images.append(striped)

        return images

def normalizetocc(point, maxvalues):
    return ((point[0] * 0.1) / float(maxvalues[0]),
            (point[1] * 0.1) / float(maxvalues[1]))

def getoutput(arn, **kwargs):
    #index = self.arn.numtf + self.arn.numrec
    index = 0
              #print leftlim
    output = arn.effectorhist[index][-1]
              #if( arn.effectorhist[index][-1] -
              #   arn.effectorhist[index][-1001]) <= 0:
              #      output = 1
              #else:
              #       output = 0
    return int(output * 255)


if __name__ == '__main__':
        arnconfigfile = 'configfiles/arnsim.cfg'
        log.setLevel(logging.DEBUG)
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(open(arnconfigfile))

        proteins=[]
        nump = 0
        try:
            f = open(sys.argv[1], 'r')
            genome = BitStream(bin=f.readline())
            arnet = ARNetwork(genome, cfg)
        except:
            while nump < 4 or nump > 32 or numeff == 0 or not arnet.receptors:
                genome = BitStream(float=random.random(), length=32)
                for i in range(cfg.getint('default','initdm')):
                    genome = dm_event(genome,
                                      .02)

                    arnet = ARNetwork(genome, cfg, numinputs = 2,numoutputs=1)
                    nump = len(arnet.promlist)
                    numeff = len(arnet.effectors)

        offspring = None
        pop = [Cell(cfg,genome,numinputs=2,numoutputs=1 )]*4
        #eval_ = bindparams(cfg, evaluatecircuit)
        plot_ = bindparams(cfg, plotindividual)

        #pop = zip()(themother, eval_(themother)),
         #      (None,0)]

        theApp = gui.App()
        theApp.images = render_images(pop, theApp)
        theApp.on_execute()
        selected = []
        while theApp._running:
            if theApp.selected:
                selected.append(pop[theApp.selected[-1]])
                plot_(selected[-1].genotype)
            pop = []
            while len(pop) < 4:
                offspring = Cell(cfg,
                                 bitflipmutation(selected[-1].genotype.code,
                                                 .01),
                                 numinputs=2,numoutputs=1)
                while (not offspring.phenotype.effectors or
                       not offspring.phenotype.receptors):
                    offspring = Cell(cfg, bitflipmutation(
                            selected[-1].genotype,.01),numinputs=2,numoutputs=1)
                pop.append(offspring)
            theApp.images = render_images(pop, theApp)
            theApp.unpause()

                #for p in arnet.proteins: print p


        #f = open('genome.save','w')
        #f.write(selected[-1].genotype.code.bin)
        #f.close
        #print genome.bin
        plot_(selected[-1].genotype)
