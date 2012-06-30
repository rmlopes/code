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

numclasses = 3
data =[map(lambda t: float(t), l.split(',')) 
       for l in open('datafiles/iris.txt').readlines()]
print data[0]

everything = [(d[4],d[:4]) for d in data]
defaultdiscard = (0,0,1,1)

def evaluate(circuit, test = False):
    if len(circuit) < 4:
        return defaultdiscard
    mainmod = __import__('__main__')
    #workingset = globals()['testset'] if test else globals()['trainset']
    workingset = getattr(mainmod,
                         'testset') if test else getattr(mainmod,'trainset')
    tp = 0 
    fp = 0 
    tn = 0 
    fn = 0
    for c,feats in workingset:
        try:
            result = evaluatecircuit(circuit, nnlikefun,
                                     dict(), *feats)
        except:
            return defaultdiscard
        
        result = 1 if result > 0 else 0
        #if test:
        #print c, result
            
        if c and result :
            tp += 1
        elif c and not result:
            fn += 1
        elif not c and not result:
            tn += 1
        elif not c and result:
            fp += 1

	    #print sample(results,10)
    return tp, tn, fp, fn    

def wrapevaluate(circuit, test = False):
    results = evaluate(circuit, test)
    if test:
        print "%i\t%i\t%i\t%i" % results
    return 1 - mcc(*results)

if __name__ == '__main__':
    import sys
    p  = ClassifProb(wrapevaluate,len(everything[0][1]))
    edw = EvoDevoWorkbench(sys.argv[1],p,buildcircuit,ReNCoDeAgent)
    
    cl = int(sys.argv[2])
    everything = map(lambda x: (1,x[1]) if x[0] == cl else (0,x[1]),
                 everything)

    random.shuffle(everything)
    splitindx = int(len(everything)*.7)
    trainset = everything[:splitindx]
    testset = everything[splitindx:]

    edw.run()
    testresult = wrapevaluate(edw.best.phenotype,True)
    print testresult

