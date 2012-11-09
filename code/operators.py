from bitstring import BitStream
from utils import *
import random

def transposon(code, *args):
    tsize = int(args[0])
    copypos = random.randint(1, len(code)-tsize)
    insertpos = random.randint(1, len(code)-tsize)
    code.insert(code[copypos:copypos+tsize],insertpos)
    return code

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

### for now just one point
def npoint_xover(code1,code2,*args):
    maxl = min(len(code1),len(code2))
    p = int(random.random() * maxl)
    o1 = BitStream(bin=(code1[:p]+code2[p:]).bin)
    o2 = BitStream(bin=(code2[:p]+code1[p:]).bin)
    return o1, o2

if __name__ == '__main__':
    a = BitStream('0b' + '1'*20)
    print a.bin
    transposon(a, 10)
    print a.bin
    junk(a, 5)
    print a.bin
    delete(a, 15)
    print a.bin

    a = BitStream('0b' + '1'*10)
    b = BitStream('0b' + '0'*20)
    c,d = npoint_xover(a,b)
    print c.bin
    print d.bin