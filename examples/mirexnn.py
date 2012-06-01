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
workingset = sample(zipped,int(0.7*len(zipped)))

class MIREX_NN(ReNCoDeProb):
    extrafuns = ['exp_','log_','tanh_','sin_','cos_','sinh_','cosh_','tan_']
    def __init__(self, evaluate):
        self.labels = None
        ReNCoDeProb.__init__(self,evaluate)
        self.funs.extend(self.extrafuns)
        self.terms.extend(['inputs[1]','inputs[2]'])
        self.arity.update(zip(self.extrafuns,[1]*len(self.extrafuns)))        
    
def nnlikefun(mapped, node_inputs, inputs):
	if not node_inputs:
            return eval(mapped)
	mainmod = __import__('__main__')
	if len(node_inputs) == 1:
            return getattr(mainmod, mapped)(node_inputs[0])
        if mapped in MIREX_NN.extrafuns:
            return getattr(mainmod, mapped)(node_inputs)
        return reduce(lambda m,n: getattr(mainmod, mapped)(m,n),
                      node_inputs)

def nn(inputs, weights1, weights2):
    hidenum = len(weights2)
    outnum = len(weights2[0])
    Z,Y = [],[]
    for i in range(hidenum):
        Z.append(sigmoid(inputs, weights1[:,i]))
    for i in range(outnum):
        Y.append(sigmoid(Z,weights2[:,i]))
    return Z,Y


def sigmoid(inputs,weights):
    return 1.0 / (1.0 + exp(-numpy.dot(inputs, weights)))


def evaluate(circuit, topology = (16,10,5), test = False):
    if len(circuit) < 4:
        return 1e4
    #numinputs = [not c[2] for c in circuit[-5:]].count(1)
    #numinputs = isinputlist.count(1)
    #if numinputs < 5:
     #   return 1e4
    #print topology
    w1 = nparray([[evaluatecircuit(circuit, nnlikefun,dict(), *[i,j,0]) 
                   for j in range(topology[1])]
                  for i in range(topology[0])])
    def thekey(e):
        if e < -1:
            return -1
        elif e > 1:
            return 1
        else:
            return e

    #print w1[1,:10]
    #w1 = nparray([map(thekey, lst) for lst in w1])
    for i in range(len(w1)):
        for j in range(len(w1[0])):
            if w1[i][j] < -1 or w1[i][j] > 1:
                return 1e4
    #print w1[1,:10]

    w2 = nparray([[evaluatecircuit(circuit, nnlikefun,dict(), *[i,j,1]) 
                   for j in range(topology[2])]
                  for i in range(topology[1])])
    #w2 = nparray([map(thekey, lst) for lst in w2])
    for i in range(len(w2)):
        for j in range(len(w2[0])):
            if w2[i][j] < -1 or w2[i][j] > 1:
                return 1e4
    
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
        
        resultv = nn(feats,w1,w2)
        resultz = zip(range(1,6),resultv[1])
        #print resultz
        rclass = max(resultz, key=lambda x: x[1])[0]
        
        #if result < 0 or result > 1.0 or isnan(result):
         #  return 1e4
        #rclass = int(result*(numclasses-1))+1
        #rclass = int(result % (numclasses-1)) + 1
        #results.append(rclass)
        if rclass == c:
            ok += 1
       
    #print sample(results,10)
    return  len(workingset) - ok    

if __name__ == '__main__':
    import sys
    p  = MIREX_NN(evaluate)
    edw = EvoDevoWorkbench(sys.argv[1],p,buildcircuit,ReNCoDeAgent)
    edw.run()
    
    testresult = evaluate(edw.best.circuit,test=True)
    print 'Generalization: '
    print testresult

