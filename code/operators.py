import random
from functools import partial
from bitstring import BitStream
#from utils import *
import logging

log = logging.getLogger(__name__)

def transposon(code, *args, **kwargs):
    tsize = int(args[0])
    copypos = random.randint(1, len(code)-tsize)
    transp = code[copypos:copypos+tsize]
    try:
        if kwargs['excise']:
            del code[copypos:copypos+tsize]
    except KeyError: pass
    insertpos = random.randint(1, len(code)-tsize)
    code.insert(transp, insertpos)
    return code

movingtransposon = partial(transposon,
                           excise = True)

def junk(code, *args):
    tsize = int(args[0])
    instream = BitStream(bin='0b'+'0'*tsize)
    insertpos = random.randint(1, len(code)-tsize)
    code.insert(instream,insertpos)
    return code

def delete(code, *args):
    tsize = int(args[0])
    delpos = random.randint(1, len(code)-tsize)
    del code[delpos:delpos+tsize]
    return code

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

def genecut_xover(p1,p2, *args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    #maxl = min(len(code1),len(code2))
    cutpts1 = random.sample(p1.genotype.promlist,2)
    cutpts1.sort()
    log.debug(str(cutpts1))
    cutpts2 = random.sample(p2.genotype.promlist,2)
    cutpts2.sort()
    log.debug(str(cutpts2))
    #p = [int(random.random() * maxl)]
    #p.append(int(random.random() * maxl))
    #p.sort()
    #p1,p2 = p
    excite_offset = p1.genotype.excite_offset
    cutpts1 = map(lambda x: x - excite_offset, cutpts1)
    cutpts2 = map(lambda x: x - excite_offset, cutpts2)
    o1 = BitStream(bin=(code1[:cutpts1[0]]+
                        code2[cutpts2[0]:cutpts2[1]]+
                        code1[cutpts1[1]:]).bin)
    o2 = BitStream(bin=(code2[:cutpts2[0]]+
                        code1[cutpts1[0]:cutpts1[1]]+
                        code2[cutpts2[1]:]).bin)
    return o1, o2

def onepoint_xover(p1,p2,*args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    maxl = min(len(code1),len(code2))
    p = int(random.random() * maxl)
    o1 = BitStream(bin=(code1[:p]+code2[p:]).bin)
    o2 = BitStream(bin=(code2[:p]+code1[p:]).bin)
    return o1, o2

def uniform_xover(p1, p2, *args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    l = min(len(code1),len(code2))
    mask = reduce(lambda x,y: "%s%s"%(x,y),
                 [int(round(random.random()))
                  for i in range(l)])
    bitmask = BitStream(bin=mask)
    #print bitmask.bin
    o1,o2 = '',''
    for i in range(len(bitmask.bin)):
        if bitmask[i]:
            o1 += code1.bin[i]
            o2 += code2.bin[i]
        else:
            o1 += code2.bin[i]
            o2 += code1.bin[i]
    return (BitStream(bin = o1), BitStream(bin = o2))

#FIXME: this operator results in a drastic reduction of genome size
# leading to bad evolution
def unigenecut_xover(p1, p2, *args):
    code1 = p1.genotype.code
    code2 = p2.genotype.code
    l = min(len(p1.genotype.promlist),len(p2.genotype.promlist))
    mask = reduce(lambda x,y: "%s%s"%(x,y),
                  [int(round(random.random()))
                   for i in range(l)])
    bitmask = BitStream(bin=mask)
    #print bitmask.bin
    o1,o2 = '',''
    for i in range(len(bitmask.bin)):
        if bitmask[i]:
            o1 += code1.bin[p1.genotype.promlist[i]-88:
                            p1.genotype.promlist[i]+169]
            o2 += code2.bin[p2.genotype.promlist[i]-88:
                            p2.genotype.promlist[i]+169]
        else:
            o1 += code2.bin[p2.genotype.promlist[i]-88:
                            p2.genotype.promlist[i]+169]
            o2 += code1.bin[p1.genotype.promlist[i]-88:
                            p1.genotype.promlist[i]+169]
    return (BitStream(bin = o1), BitStream(bin = o2))

if __name__ == '__main__':
    log.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    log.addHandler(ch)
    #test inner variation operators
    a = BitStream('0b' + '00011100011100011100')
    print a.bin
    transposon(a, 5)
    print a.bin
    junk(a, 5)
    print a.bin
    delete(a, 15)
    print a.bin

    #test crossover operators
    class Agent: pass
    a = Agent()
    a.genotype = Agent()
    a.genotype.code = BitStream('0b' + '1'*10)
    #a =
    b = Agent()
    b.genotype = Agent()
    b.genotype.code = BitStream('0b' + '0'*20)
    #b =
    c,d = twopoint_xover(a,b)
    print c.bin
    print d.bin
    e,f = uniform_xover(a,b)
    print e.bin
    print f.bin

    #test gene oriented crossover
    import ConfigParser
    from rencode import DMAgent, ReNCoDeProb
    arncfg = ConfigParser.ConfigParser()
    arncfg.readfp(open('../configfiles/arnsim-rencode.cfg'))
    p = ReNCoDeProb(None)
    a = DMAgent(arncfg, p)
    b = DMAgent(arncfg, p)
    print "A: ", a.genotype.promlist
    print "B: ", b.genotype.promlist
    c,d = map(lambda o: DMAgent(arncfg,p,o),genecut_xover(a,b))
    print "C: ", c.genotype.promlist
    print "D: ", d.genotype.promlist
    e,f = map(lambda o: DMAgent(arncfg,p,o),unigenecut_xover(a,b))
    print "E: ", e.genotype.promlist
    print "F: ", f.genotype.promlist