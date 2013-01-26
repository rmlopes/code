#from code import operators
from code.evodevo import *
from code.operators import *
from math import *
from code.utils.mathlogic import *
from code.rencode2 import *
from code.rencode import evaluatecircuit, ReNCoDeProb, regressionfun
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

def evaluate(phenotype, test = False, relax = False):
    circuit = phenotype.getcircuit()
    if len(circuit) < 4 and not relax:
        return 1.0

    mainmod = __import__('__main__')
    #workingset = globals()['testset'] if test else globals()['trainset']
    workingset = getattr(mainmod,
                         'testset') if test else getattr(mainmod,'trainset')

    confusion = nparray([[0] * 5
                         for i in range(5)])
    for c,feats in workingset:
        #try:
        results = evaluatecircuit(circuit, regressionfun,
                                  dict(), *feats, nout=5)

        to_order = zip(results,range(1,6))
        to_order.sort(key=lambda x: x[0], reverse = True)
        predicted = to_order[0][1]
        confusion[c-1][predicted-1] += 1

    partials = confusion_f1(confusion)
    #print confusion
    partials = filter(lambda x: not isnan(x), partials)
    #print partials
    if test:
        return confusion
    return 1 - sum(partials)/5.0

def confusion_f1(confusion):
    totalf1 = tuple()
    for i in range(confusion.shape[0]):
        p = confusion_precision(confusion,i)
        r = confusion_recall(confusion,i)
        totalf1 += (2*p*r/(p+r),)
    return totalf1

def confusion_precision(confusion, target):
    '''target is 0 indexed'''
    numclasses = confusion.shape[0]
    tp = confusion[target][target]
    e = sum([c for c,i in zip(confusion[:,target],range(numclasses))
             if i != target])
    return tp/float(tp+e)


def confusion_recall(confusion, target):
    '''target is 0 indexed'''
    numclasses = confusion.shape[0]
    tp = confusion[target][target]
    e = sum([c for c,i in zip(confusion[target,:],range(numclasses))
             if i != target])
    return tp/float(tp+e)

globalconfusion = nparray([[0]*numclasses
                           for i in range(numclasses)])

def runfold(i, datasets, edw, globalconfusion):
    print "Running fold %i"%i
    mainmod = __import__('__main__')
    testset = []
    trainset = []
    for j in range(len(datasets)):
        if j != i:
            trainset.extend(datasets[j])
        if j == i:
            testset.extend(datasets[j])
    setattr(mainmod,'trainset',trainset)
    setattr(mainmod,'testset',testset)
    edw.run()
    confusion = evaluate(edw.best.phenotype, True)
    globalconfusion += confusion
    return confusion_f1(confusion)

class CProb(ReNCoDeProb):
    extrafuns = ['add_','sub_','mul_','div_']
    def __init__(self, evaluate, numinputs):
        self.nout = 5
        self.ninp = numinputs
        self.labels = None
        ReNCoDeProb.__init__(self,evaluate,printf=printmultiplecircuit)
        self.terms.extend(["inputs[%i]"%i
                           for i in range(1,numinputs)])
        self.funs.extend(self.extrafuns)
        self.arity.update(zip(self.extrafuns,[0]*len(self.extrafuns)))

def createfolds(data, numfolds = 3):
    #zipped1 = map(lambda x: (1,x[1],(x[0],x[2])) if x[0] == target else (0,x[1],(x[0],x[2])),
     #             data)
    #cset = filter(lambda x: x[0] == target, zipped1)
    #otherset = filter(lambda x: x[0] != target, zipped1)
    #workset = cset + random.sample(otherset, len(cset))
    random.shuffle(data)
    foldsize = int(len(data)/numfolds)
    folds = []
    for i in range(numfolds-1):
        start = i*foldsize
        folds.append(data[start:start+foldsize])
    folds.append(data[(numfolds-1)*foldsize:])
    return folds

if __name__ == '__main__':
    import sys
    try:
        if sys.argv[2]  == 'pca': #nopca
            zipped = zip(allclasses, projected)#, range(len(allclasses)))
        elif sys.argv[2] == 'original':
            zipped = zip(allclasses, allfeatures)#, range(len(allclasses)))

        numfolds = int(sys.argv[3])
        print "NumFeatures: %i "% (len(zipped[0][1]))

        #fit = sys.argv[4]
    except:
        print "usage: python -m examples.mirex <config> <feat:pca|original> <numfolds>"
        print "example: python -m examples.mirex test.cfg pca 3"
        exit(0)

    datasets = []
    datasets = createfolds(zipped, numfolds)
    #for i in range(1,6):
     #   datasets.append(createdataset(zipped,i))

    problem  = CProb(evaluate,len(zipped[0][1]))

    #problem.funs = ['add_','sub_','mul_','div_','min_','max_']
    edw = EvoDevoWorkbench(sys.argv[1],problem,RndAgent)

    vresults = 0
    for i in range(numfolds):
        v = runfold(i, datasets, edw, globalconfusion)
        avgf1 = sum(v)/5.0
        log.critical("%f\t%f\t%f\t%f\t%f\t%f"% (tuple(v)+(avgf1,)))
        vresults += avgf1
    ##v = runfolds(edw)

    print 'CONFUSION:'
    for i in range(5):
        print "%i\t%i\t%i\t%i\t%i" % tuple(globalconfusion[i,:])

    log.critical("avg. f1 of the folds: %f", vresults/numfolds)
