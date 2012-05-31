#from code import operators
from code.evodevo import *
from code.operators import *
from math import *
from code.utils.mathlogic import *
from code.rencode import *
from random import sample
#import matplotlib.mlab as mlab
from numpy import array as nparray
import numpy

numclasses = 5
allclasses = nparray([int(c) 
                      for c in open('datafiles/MIREXclasses.txt').readlines()])
allfeatures =nparray([map(lambda t: float(t), l.split(',')) 
                      for l in open('datafiles/FMnorm.txt').readlines()])

#pcafeats = mlab.PCA(allfeatures)
#projected = pcafeats.Y[:,:16]#pcafeats.project(allfeatures)
projected = numpy.load('datafiles/projectedfeat.npy')
myitercount = 0
zipped = zip(allclasses, projected)
workingset = sample(zipped,int(0.3*len(zipped)))

class MIREX(ReNCoDeProb):
    extrafuns = ['exp_','log_','tanh_','sqrt_']
    def __init__(self, evaluate):
        self.labels = None
        ReNCoDeProb.__init__(self,evaluate)
        self.terms.extend(["inputs[%i]"%i 
                           for i in range(1,len(zipped[0][1]))])
        self.funs.extend(self.extrafuns)
        self.arity.update(zip(self.extrafuns,[1]*len(self.extrafuns)))
    
def nnlikefun(mapped, node_inputs, inputs):
	if not node_inputs:
            return eval(mapped)
	mainmod = __import__('__main__')
	if len(node_inputs) == 1:
            return getattr(mainmod, mapped)(node_inputs[0])
        if mapped in MIREX.extrafuns:
            return getattr(mainmod, mapped)(node_inputs)
        return reduce(lambda m,n: getattr(mainmod, mapped)(m,n),
                      node_inputs)

def evaluate(circuit, test = False):
    if len(circuit) < 4:
        return 1e4
    #numinputs = [not c[2] for c in circuit[-5:]].count(1)
    #numinputs = isinputlist.count(1)
    #if numinputs < 5:
     #   return 1e4
    ok = 0
    last1 = 0
    last2 = 0
    if test :
        workingset = zipped
    else:
        #if  globals()['myitercount'] < edw.itercount:
        #    globals()['myitercount'] = edw.itercount
        #    globals()['workingset'] = sample(zipped,int(0.3*len(zipped)))
        workingset = globals()['workingset']
    
    results = []
    for c,feats in workingset:
        try:
            result = evaluatecircuit(circuit, nnlikefun,
                                     dict(), *feats)
        except:
            return len(workingset)
        if result < 0 or result > 1.0 or isnan(result):
           return 1e4
        rclass = int(result*(numclasses-1))+1
        #rclass = int(result % (numclasses-1)) + 1
        #results.append(rclass)
        if rclass == c:
            ok += 1
       
    #print sample(results,10)
    return  len(workingset) - ok    

if __name__ == '__main__':
    import sys
    total = len(allclasses) 
    splitindex = int(total*.3)
    #classes,features = zip(*(sample(zipped, splitindex)))
    print 'Evolving with %i examples' % splitindex

    p  = MIREX(evaluate)
    edw = EvoDevoWorkbench(sys.argv[1],p,buildcircuit,ReNCoDeAgent)
    edw.run()
    
    testresult = evaluate(edw.best.circuit,True)
    print 'Generalization: '
    print testresult

