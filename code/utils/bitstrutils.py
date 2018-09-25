from bitstring import BitArray
import random
import numpy
from math import floor,ceil
from bitarray import *

def dm_event(code, mutation_rate):
    code += bitflipmutation(bitarray(code),mutation_rate)
    return code

def XORmatchingDual(pcode,ecode,inhcode,threshold):
    ematch = pcode ^ ecode
    inhmatch = pcode ^ inhcode
    result = []
    if ematch.count(1) < threshold: result.append(0)
    else : result.append(ematch.count(1))

    if inhmatch.count(1) < threshold: result.append(0)
    else : result.append( inhmatch.count(1) )
    return result

def XORmatching(pcode,ecode,threshold):
    match = (pcode ^ ecode).count()
    if match < threshold: return 0
    return match

def bitflipmutation(code, mutrate = 0.01, **kwargs):
    rnds = [1 if random.random() <= mutrate else 0 for i in xrange(code.length())]
    r = bitarray(rnds)
    return code ^ r

def epibitflipmutation(code, mask, rate):
    for i in range(0,len(code)):
        if random.random() < rate and not mask[i]:
            code.invert(i)
    return code

def applymajority__(code, chunk_size):
    chunknum = len(code) / chunk_size
    lim = chunknum // 2 #floor is too slow
    bits=''
    for i in range(0,chunk_size):
        codeset = [code[i + j*chunk_size]
                   for j in range(0,chunknum)]
        c = codeset.count()
        bits += '1' if c > lim else '0'
    return bits

def applymajority(code, chunk_size):
    lim = ceil(code.length() / float(chunk_size * 2))
    l = ['1' if code[i::chunk_size].count()>=lim else '0' for i in xrange(chunk_size)]
    return ''.join(l)

def id_generator(size=32, chars=['0','1']):
    return ''.join(random.choice(chars) for x in range(size))

def getrndstr(size):
    return BitArray(int=random.randint(0,2**(size-1)), length=size).bin



if __name__=='__main__':
    import cProfile, pstats, StringIO
    from bitarray import *
    pr = cProfile.Profile()
    pr.enable()
    random.seed(1234)
    c1 = bitarray(getrndstr(128))
    print c1.to01()
    c2 = bitarray(getrndstr(128))
    for i in range(10000):
        #m = XORmatching_(c1,c2,0)
        #m = applymajority(c1,32)
        m = bitflipmutation(c1,0.01)
        if i == 0:
            print m.to01()
    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()
    exit(0)
    i=0
    while i<1:
        ps=[]
        mutbs=[]
        for b in bs:
            mutbs.append(bitflipmutation(BitStream(b),0.5))
        for mb in mutbs:
            signature = applymajority(mb,32)
            ps.append( BitStream(bin=signature))
        i+=1

        print XORmatching(BitStream('0b110'),BitStream('0b111'),0)
        print XORmatchingDual(BitStream('0b000'),
                              BitStream('0b100'),BitStream('0b111'),0)
        print bs[0].bin[0:32]
        print mutbs[0].bin[0:32]
        print ps[0].bin
