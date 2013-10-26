from bitstring import BitStream
import random
from math import floor

def dm_event(code, mutation_rate):
    temp = BitStream(bin=code.bin)
    return code + bitflipmutation(temp, mutation_rate)

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
    ematch = pcode ^ ecode
    match = ematch.count(1)
    if match < threshold: return 0
    return match

def bitflipmutation(code, mutrate = 0.01, **kwargs):
    for i in range(0,len(code)):
        if random.random() < mutrate:
            code.invert(i)
    return code

def epibitflipmutation(code, mask, rate):
    for i in range(0,len(code)):
        if random.random() < rate and not mask[i]:
            code.invert(i)
    return code

def applymajority(code, chunk_size):
    chunknum = len(code) / chunk_size
    lim = floor(chunknum/2)
    bits=''
    for i in range(0,chunk_size):
        codeset = [code[i + j*chunk_size]
                   for j in range(0,chunknum)]
        c = codeset.count(1)
        bits += '1' if c > lim else '0'
    return bits

if __name__=='__main__':
    bs = [BitStream(bin='00100000110000100000100000111000'*5)]
    #[BitStream(bin='011101011110001110001110001101110000011111000011'*5)]

    i=0
    while i<1:#True:
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
