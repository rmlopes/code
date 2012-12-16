import ConfigParser, os
import logging
import random
import numpy as np
from numpy import array as nparray
from functools import partial
from bitstring import *
from math import exp
import sys
from sys import stdout
from subprocess import call
from code.utils import *
from code.utils.bitstrutils import *
#from code.extendedarn import *
from code.cellcode import *
import cPickle as pickle
#from code.arnmiguel import *


log = logging.getLogger(__name__)


if __name__ == '__main__':
        arnconfigfile = sys.argv[2]
        #'configfiles/arnsim-miguel.cfg'
        log.setLevel(logging.DEBUG)
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(open(arnconfigfile))
        evalfun = evaluatecircuit
        #mapfun=getoutputp0p1)
        p = CellProb(evalfun, 3, 1)
        #try:
            #genome = BitStream(bin=f.readline())
            #arnet = ARNetwork(genome, cfg, problem = p)
        c = pickle.load(open(sys.argv[1], 'r'))
        c.reset()
            #print arnet.proteins
        #except:
         #   print 'Failed to load test individual.'

        p.eval_ = bindparams(cfg, evaluatecircuit)
        p.plot_ = bindparams(cfg, plotindividual)
        #c = Cell(cfg,genome, problem = p )
        for i in range(2):
                fit = p.eval_(c.phenotype, True, shuffle = False)
                print 'FIT: ', fit
        p.plot_(c.phenotype)
        for i in range(6):
                fit = p.eval_(c.phenotype, True, shuffle = True)
                print 'FIT: ', fit
        while True:
                fit = p.eval_(c.phenotype, True, shuffle = False)
                print 'FIT: ', fit
        #fit = eval_(c, shuffle = False)
        p.plot_(c.phenotype)
        print 'FIT: ', fit
