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
from code.extendedarn import *
from code.cellcode import *

log = logging.getLogger(__name__)


if __name__ == '__main__':
        arnconfigfile = 'configfiles/arnsim.cfg'
        log.setLevel(logging.DEBUG)
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(open(arnconfigfile))
        proteins=[]
        nump = 0
        prob_inp=[.0,.0,.0]
        p = CellProb(evaluatecircuit, 3, 1)
        try:
            f = open(sys.argv[1], 'r')
            genome = BitStream(bin=f.readline())
            arnet = ARNetwork(genome, cfg, problem = p)
        except:
            while nump < 4 or nump > 32 or numeff == 0 or not arnet.receptors:
                genome = BitStream(float=random.random(), length=32)
                for i in range(cfg.getint('default','initdm')):
                    genome = dm_event(genome,
                                      .02)

                    arnet = ARNetwork(genome, cfg,problem = p)
                    nump = len(arnet.promlist)
                    numeff = len(arnet.effectors)

        eval_ = bindparams(cfg, evaluatecircuit)
        plot_ = bindparams(cfg, plotindividual)
        c = Cell(cfg,genome, problem = p )
        i = 0
        s = False
        while True:
                fit = eval_(c, True, shuffle = False)
                print 'FIT: ', fit
                i += 1
        #fit = eval_(c, shuffle = False)
        plot_(c.phenotype)
        print 'FIT: ', fit
