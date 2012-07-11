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

projected = numpy.load('datafiles/projectedfeat-01.npy')

def multiplexadf(inputs, circuit):
    return evaluatecircuit(circuit, nnlikefun,
                           dict(), *inputs)

from mirexadf import createADFs, runevo, runevofolds, fmeas_eval
from iris import wrapevaluate

if __name__ == '__main__':
    import sys
    try:
        if sys.argv[3]  == 'pca': #nopca
            zipped = zip(allclasses, projected)
        elif sys.argv[3] == 'original':
            zipped = zip(allclasses, allfeatures)

        numfolds = int(sys.argv[4])
        print "NumFeatures: %i "% (len(zipped[0][1]))

    #Class by command line
        c = int(sys.argv[2])
    except:
        print "usage: python -m examples.mirex <config> <classnum> <features> <numfolds>"
        print "example: python -m examples.mirex test.cfg 1 pca 3"
        exit(0)

    #create evolutionary workbench
    problem  = ClassifProb(wrapevaluate,len(zipped[0][1]))
    #problem.terms = ['adf0']
    edw = EvoDevoWorkbench(sys.argv[1],problem,buildcircuit,ReNCoDeAgent)

    #create class dataset
    zipped1 = map(lambda x: (1,x[1]) if x[0] == c else (0,x[1]),
                  zipped)

    cset = filter(lambda x: x[0] == 1, zipped1)
    otherset = filter(lambda x: x[0] == 0, zipped1)
    #used when no folds
    trainset = cset[:int(.7*len(cset))]
    trainset.extend(random.sample(otherset, len(trainset)))
    random.shuffle(trainset)
    testset = cset[int(.7*len(cset)):]
    testset.extend(random.sample(otherset, len(testset)))
    random.shuffle(testset)
    workset = cset + random.sample(otherset, len(cset))
    print "workset size: %i" % (len(workset),)

    #no folds
    if(numfolds == 1):
       #createADFs(edw, problem, trainset)
       edw.run()
       fit = fmeas_eval(edw.best.phenotype,True)
       print "F:\n",fit
       exit(0)
    else:
        #create folds
        random.shuffle(workset)
        #createADFs(edw, problem, workset)
        foldsize = int(len(workset)/numfolds)
        folds = []
        for i in range(numfolds-1):
            start = i*foldsize
            folds.append(workset[start:start+foldsize])
        folds.append(workset[(numfolds-1)*foldsize:])
        print "Created %i folds" % (numfolds,)

        fmeas, bestofbest = runevofolds(edw, problem, folds)



'''
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
'''
