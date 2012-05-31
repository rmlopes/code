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

if __name__ == '__main__':
    a = BitStream('0b' + '1'*20)
    print a.bin
    transposon(a, 10, 0.5, 20)
    print a.bin
    junk(a, 5, 0.5, 20)
    print a.bin
    delete(a, 15, 0.5, 20)
    print a.bin
