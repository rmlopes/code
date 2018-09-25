# test single pole performance
import single_pole
from bitstring import *
from math import exp
import sys
from sys import stdout
from subprocess import call
from code.utils import *
from code.utils.bitstrutils import *
from code.arnmiguel import *
from code.cellcode import *
from random import randint
import cPickle as pickle


TWELVE_DEGREES = 0.2094384 #radians

def factorstoinputs( factors ):
    return [(factors[0] * 4.8) - 2.4,
            (factors[1] * 2.0) - 1.0,
            (factors[2] * (TWELVE_DEGREES*2)) - TWELVE_DEGREES,
            (factors[3] * 3.0) - 1.5]

if __name__ == '__main__':
    arnconfigfile = sys.argv[1]
    #'configfiles/arnsim-miguel.cfg'
    log.setLevel(logging.DEBUG)
    cfg = ConfigParser.ConfigParser()
    cfg.readfp(open(arnconfigfile))
    proteins=[]
    nump = 0
    prob_inp=[.0]*4
    p = CellProb(single_pole.evaluate_individual, 4, 1)
    p.eval_ = bindparams(cfg, p.eval_)

    f = open(sys.argv[2]+os.getenv('SGE_TASK_ID')+'.save', 'r')

    (binstr,outidx) = pickle.load(f)
    genome = BitStream(bin=binstr)
    cell = Cell(cfg, genome, problem = p)
    cell.phenotype.output_idx = outidx
    testportions = [.05, .275, .5, .725, .95]
    import itertools
    perms = itertools.product(testportions,repeat = 4)
    inputs = map(factorstoinputs, perms)
    assert len(inputs) == 625
    results = []
    original = copy.deepcopy(cell.phenotype.ccs)
    for i in inputs:
        fit = p.eval_(cell.phenotype, i)
        if fit < 0.000000001:
            results.append(0)
        else:
            results.append(1)
        print i, ' ', fit
        cell.phenotype.reset(original)
     #print "\nInitial conditions:"
     #print "%2.4f   %2.4f   %2.4f   %2.4f" %tuple(init)

    rstr = ''
    for i in range(len(results)):
        rstr = rstr + str(results[i]) + ' '
        if (i+1) % pow(5,2) == 0:
            rstr = rstr + '\n'
    print rstr

    print "Generalization Fitness: %f" % (625 - sum(results),)
