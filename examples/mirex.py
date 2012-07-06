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
projected = numpy.load('datafiles/projectedfeat-02.npy')
    
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
    mainmod = __import__('__main__')
    terminals = set([n[1] for n in circuit if not n[2]])
    p =  getattr(mainmod,'problem')
    #print terminals
    #print p.terms
    #penalty = len(p.terms) - len(terminals)
    penalty = 0
    results = evaluate(circuit, test)
    tp, tn, fp, fn = results
    if test:
        log.critical("%i\t%i\t%i\t%i", *results)
    try: 
        f = fmeasure(tp, fp, fn)
    except ZeroDivisionError:
        f = -1

    return (1 - f)*(1 + penalty) 
    

if __name__ == '__main__':
    import sys    
    if sys.argv[3]  == 'pca': #nopca
        zipped = zip(allclasses, projected)
    elif sys.argv[3] == 'original':
        zipped = zip(allclasses, allfeatures)


    print "NumFeatures: %i "% (len(zipped[0][1]))

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
    for i in range(9):
        start = i*foldsize
        folds.append(workset[start:start+foldsize]) 
    folds.append(workset[9*foldsize:])
    print "Created %i folds" % (len(folds),)

    problem  = ClassifProb(fmeas_eval,len(zipped[0][1]))
    edw = EvoDevoWorkbench(sys.argv[1],problem,buildcircuit,ReNCoDeAgent)
        
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

