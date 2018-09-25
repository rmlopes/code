import sys
from code.utils.config import parsecmd, loadconfig
from code.evodevo import *
from code.operators import *
from code.rencode2 import *
from code.rencode import ReNCoDeProb, evaluatecircuit, mergefun
from code.utils.mathlogic import *
from code.extendedarn import bindparams, displayARNresults
import logging

log = logging.getLogger(__name__)

class BooleanProb(ReNCoDeProb):
        labels = {'and_':'AND',
                  'or_':'OR',
                  'nand':'NAND',
                  'nor':'NOR',
                  'inputs[0]':'INP'}
        funs = [ 'and_', 'or_', 'nand', 'nor','and_','or_','nand','nor']
        feedback = False
        def __init__(self, nbits, evalf, **kwargs):
                self.ninp = nbits*2
                self.nout = nbits*2
                self.terms = ['inputs[%i]'%(idx,) for idx in range(nbits*2)]
                ReNCoDeProb.__init__(self,evalf, **kwargs)




if __name__ == '__main__':
    import sys
    import os
    import cPickle as pickle
    #from pydot import graph_from_dot_data

    fname = '../../Documents/thesis/phdsupport/data/multiple/analysis/multiplier.save'
    tfile = open(fname,'r')
    #lines= tfile.readlines()
    #outfile = sys.argv[2]
    #print tfile.read()
    gstr = tfile.readlines()[0][3:-2]
    #g = pickle.loads(tfile.readlines()[0][1:])
    print gstr
    #log.setLevel(logging.DEBUG)
    #random.seed(1234*int(os.getenv('SGE_TASK_ID')))
    #p  = BooleanProb(evaluate)
    #evalf = partial(None, nbits=2)
    #mapfun = getoutputp0p1)
    cfg = loadconfig(parsecmd())
    p = BooleanProb(2,evaluatecircuit)
    a = RndAgent(config = cfg[1],problem = p, gcode=BitStream(bin=gstr))
    a.phenotype.funskel = mergefun
    print  printmultiplecircuit(a.phenotype, labels=None, arnet = a.genotype)
    print a.phenotype(*['inputs[0]','inputs[1]','inputs[2]','inputs[3]'], nout=4)
    #edw = EvoDevoWorkbench(cfg,p)
    #p.eval_ = bindparams(edw.arnconfig, p.eval_)
    #edw.run()

    #f = open('genome.save','w')
    #f.write(edw.best.phenotype.code.bin)
    #zzf.close

    #plot_ = bindparams(edw.arnconfig, plotindividual)
    #plot_(edw.best.phenotype)
    #genresult = test(ewd.best.circuit, range(2,20))
    #print 'Generalization: '
    #print genresult
