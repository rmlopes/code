import random
import logging

log = logging.getLogger(__name__)

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


def createfolds(data, numfolds = 3):
    random.shuffle(data)
    foldsize = int(len(data)/numfolds)
    folds = []
    for i in range(numfolds-1):
        start = i*foldsize
        folds.append(data[start:start+foldsize])
    folds.append(data[(numfolds-1)*foldsize:])
    return folds


TRAINSET = 'trainset'
TESTSET = 'testset'

def setcurrent(i, folds):
    log.info("Setting the test and train sets based on fold %i."%i)
    mainmod = __import__('__main__')
    testset = []
    trainset = []
    for j in range(len(folds)):
        if j != i:
            trainset.extend(folds[j])
        if j == i:
            testset.extend(folds[j])
    setattr(mainmod,TRAINSET,trainset)
    setattr(mainmod,TESTSET,testset)
