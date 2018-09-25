import random
from functools import partial
from bitstring import BitStream
from bitarray import *
#from utils import *
import logging

log = logging.getLogger(__name__)

def transposon(code, *args, **kwargs):
    tsize = int(args[0])
    copypos = random.randint(1, len(code)-tsize)
    transp = code[copypos:copypos+tsize]
    insertpos = random.randint(1, len(code)-tsize)
    c1 = code[:insertpos]
    c1 += transp
    c1 +=code[insertpos:]
    return c1

def movingtransposon(code, *args, **kwargs):
    tsize = int(args[0])
    copypos = random.randint(1, len(code)-tsize)
    transp = code[copypos:copypos+tsize]
    insertpos = random.randint(1, len(code)-tsize)
    code[insertpos:insertpos+tsize] = transp
    return code

def junk(code, *args, **kwargs):
    tsize = int(args[0])
    instream = bitarray('0'*tsize)
    insertpos = random.randint(1, len(code)-tsize)
    c = code[:insertpos]
    c += instream
    c += code[insertpos:]
    return c

def delete(code, *args, **kwargs):
    tsize = int(args[0])
    delpos = random.randint(1, len(code)-tsize)
    del code[delpos:delpos+tsize]
    return code

def onepoint_xover(p1,p2,*args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    maxl = min(len(code1),len(code2))
    p = int(random.random() * maxl)
    o1 = BitStream(bin=(code1[:p]+code2[p:]).bin)
    o2 = BitStream(bin=(code2[:p]+code1[p:]).bin)
    return o1, o2

def twopoint_xover(p1,p2, *args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    maxl = min(len(code1),len(code2))
    p = [int(random.random() * maxl)]
    p.append(int(random.random() * maxl))
    p.sort()
    p1,p2 = p
    o1 = BitStream(bin=(code1[:p1]+code2[p1:p2]+code1[p2:]).bin)
    o2 = BitStream(bin=(code2[:p1]+code1[p1:p2]+code2[p2:]).bin)
    return o1, o2

def uniform_xover(p1, p2, *args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    l = min(len(code1),len(code2))
    s = random.randint(0, 2**l - 1)
    bitmask = BitStream(uint=s,length=l)
    rbitmask = BitStream(~bitmask, length=len(bitmask))
    t1 =  code1 & bitmask
    t2 = code2 & (rbitmask)
    o1 = t1 | t2
    t1 =  code2 & bitmask
    t2 = code1 & (rbitmask)
    o2 = t1 | t2
    return o1,o2

def onepointgene_xover(p1,p2,*args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    maxl = min(len(p1.genotype.promlist),len(p2.genotype.promlist))
    eos = p1.genotype.excite_offset
    p = int(random.random() * maxl)
    o1 = BitStream(bin=(code1[:p1.genotype.promlist[p]-eos]+code2[p2.genotype.promlist[p]-eos:]).bin)
    o2 = BitStream(bin=(code2[:p2.genotype.promlist[p]-eos]+code1[p1.genotype.promlist[p]-eos:]).bin)
    return o1, o2


def twopointgene_xover(p1,p2,*args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    maxl = min(len(p1.genotype.promlist),len(p2.genotype.promlist))
    eos = p1.genotype.excite_offset
    cutpts = random.sample(range(maxl),2)
    cutpts.sort()
    realcut1 = p1.genotype.promlist[cutpts[0]]-eos , p1.genotype.promlist[cutpts[1]]-eos
    realcut2 = p2.genotype.promlist[cutpts[0]]-eos , p2.genotype.promlist[cutpts[1]]-eos
    o1 = BitStream(bin=(code1[:realcut1[0]]+code2[realcut2[0]:realcut2[1]]+code1[realcut1[1]:]).bin)
    o2 = BitStream(bin=(code2[:realcut2[0]]+code1[realcut1[0]:realcut1[1]]+code2[realcut2[1]:]).bin)
    return o1, o2


def unigene_xover(p1, p2, *args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    l = min(len(p1.genotype.promlist),len(p2.genotype.promlist))
    s = random.randint(0, 2**l - 1)
    bitmask = BitStream(uint=s,length=l)

    o1 = code1.bin[:p1.genotype.promlist[0]-88]
    o2 = code2.bin[:p2.genotype.promlist[0]-88]
    for i in range(len(bitmask.bin)):
        if bitmask[i]:
            o1 += code1.bin[p1.genotype.promlist[i]-88:
                            p1.genotype.promlist[i]+168]
            o2 += code2.bin[p2.genotype.promlist[i]-88:
                            p2.genotype.promlist[i]+168]
        else:
            o1 += code2.bin[p2.genotype.promlist[i]-88:
                            p2.genotype.promlist[i]+168]
            o2 += code1.bin[p1.genotype.promlist[i]-88:
                            p1.genotype.promlist[i]+168]

    o1 +=  code1.bin[p1.genotype.promlist[-1]+168:]
    o2 +=  code2.bin[p2.genotype.promlist[-1]+168:]
    return (BitStream(bin = o1), BitStream(bin = o2))

def genecopy(code, *args, **kwargs):
    arnet = kwargs['arnet']
    choice = random.choice(arnet.promlist)
    #TODO: remove hard coded constants
    transp = code[choice-88:choice+168]
    #always at the end because it is indifferent
    code.append(transp)
    return code

def genedelete(code, *args, **kwargs):
    arnet = kwargs['arnet']
    choice = random.choice(arnet.promlist)
    #TODO
    del code[choice-88:choice+168]
    return code

if __name__ == '__main__':
    import ConfigParser
    from rencode import DMAgent, ReNCoDeProb
    arncfg = ConfigParser.ConfigParser()
    arncfg.readfp(open('configfiles/arnsim-rencode.cfg'))
    log.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    log.addHandler(ch)
    p = ReNCoDeProb(None)
    #test inner variation operators
    a = bitarray('00011100011100011100')
    print a.to01()
    a1 = transposon(bitarray(a), 5)
    print a1.to01()
    a2 = junk(bitarray(a), 5)
    print a2.to01()
    a3 = delete(bitarray(a), 15)
    print a3.to01()
    a4 = movingtransposon(bitarray(a), 5)
    print a4.to01()
    exit(0)

    #test genecopy and delete
    a = DMAgent(arncfg, p)
    print "Parent: ", a.genotype.promlist
    b = DMAgent(arncfg,p,genecopy(BitStream(a.genotype.code),
                                  arnet = a.genotype))
    print "Offspring: ", b.genotype.promlist
    c = DMAgent(arncfg,p,genedelete(BitStream(a.genotype.code),
                                  arnet = a.genotype))
    print "Offspring: ", c.genotype.promlist
    #exit(0)
    #test crossover operators
    class Agent: pass
    a = Agent()
    a.genotype = Agent()
    a.genotype.code = BitStream('0b' + '1'*20)
    #a =
    b = Agent()
    b.genotype = Agent()
    b.genotype.code = BitStream('0b' + '0'*20)
    #b =
    c,d = twopoint_xover(a,b)
    print c.bin
    print d.bin
    print "UNIFORM:"
    e,f = uniform_xover(a,b)
    print e.bin
    print f.bin
    #exit(0)
    #test gene oriented crossover


    a = DMAgent(arncfg, p)
    b = DMAgent(arncfg, p)
    print "A: ", a.genotype.promlist
    print "B: ", b.genotype.promlist
    c,d = map(lambda o: DMAgent(arncfg,p,o),twopointgene_xover(a,b))
    print "C: ", c.genotype.promlist
    print "D: ", d.genotype.promlist
    exit(0)
    e,f = map(lambda o: DMAgent(arncfg,p,o),unigenecut_xover(a,b))
    print "E: ", e.genotype.promlist
    print "F: ", f.genotype.promlist
