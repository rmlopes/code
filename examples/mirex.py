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
import logging

log = logging.getLogger('mirex')

numclasses = 5
allclasses = nparray([int(c) 
                      for c in open('datafiles/MIREXclasses.txt').readlines()])
allfeatures =nparray([map(lambda t: float(t), l.split(',')) 
                      for l in open('datafiles/FMnorm.txt').readlines()])

#pcafeats = mlab.PCA(allfeatures)
#projected = pcafeats.Y[:,:16]#pcafeats.project(allfeatures)
#projected = numpy.load('datafiles/projectedfeat-01.npy')
zipped = zip(allclasses, allfeatures)#projected)
    
def evaluatemulticlass(circuit, test = False):
    if len(circuit) < 4:
        return 1e4
    #numinputs = [not c[2] for c in circuit[-5:]].count(1)
    #numinputs = isinputlist.count(1)
    #if numinputs < 5:
     #   return 1e4
    ok = 0
    last1 = 0
    last2 = 0

    workingset = globals()['testset'] if test else globals()['trainset']

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

from iris import evaluate

def fmeas_eval(circuit, test = False):
    results = evaluate(circuit, test)
    tp, tn, fp, fn = results
    if test:
        log.critical("%i\t%i\t%i\t%i", results[0],results[1],results[2],results[3])
    try: 
        f = fmeasure(tp, fp, fn)
    except ZeroDivisionError:
        f = 0
    return 1 - f 
    

if __name__ == '__main__':
    import sys    
    splitindex1 = int(.7*len(zipped))
    #splitindex2 = int(.2*len(zipped)) + splitindex1
    #Class by command line
    c = int(sys.argv[2])
    zipped1 = map(lambda x: (1,x[1]) if x[0] == c else (0,x[1]),
             zipped)
    #random.shuffle(zipped1)
    
    cset = filter(lambda x: x[0] == 1, zipped1)
    otherset = filter(lambda x: x[0] == 0, zipped1)
    workset = cset + random.sample(otherset, len(cset))
    random.shuffle(workset)
    print "workset size: %i" % (len(workset),)

    foldsize = int(.1*len(workset))
    folds = []
    for i in range(10):
        start = i*foldsize
        folds.append(workset[start:start+foldsize]) 
    print "Created %i folds" % (len(folds),)

    p  = ClassifProb(fmeas_eval,len(zipped[0][1]))
    edw = EvoDevoWorkbench(sys.argv[1],p,buildcircuit,ReNCoDeAgent)
        
    vresults = 0
    for i in range(10):
        testset = folds[i]
        trainset = []
        for j in range(10):
            if j != i:
                trainset.extend(folds[j])
        
        print "Train size is %i and Test size is %i" % (len(trainset),len(testset)) 
    
        edw.run()
        vresults += fmeas_eval(edw.best.phenotype,True)
        
    log.critical("%f", vresults/10)

