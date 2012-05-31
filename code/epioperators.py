from bitstring import BitStream
from bitstringutils import *
import random

def dm_op(code, oprate, mutrate):
    if random.random() < rate:
        code = dm_event(code, mutrate)
    return code

#IMPLEMENT BETTER MORE GENERALIST HUB
#THIS IS CRAP FOR THIS MODEL
def dualophub(code, mask, op1, op2, rate1, rate2,mutfun):
    r = random.random()
    if r < rate1: 
        op1(code,epicode)
    elif r < rate1+rate2: 
        op2(code,epicode)
    else:
        mutfun(code,mask)
    return code

def transposon(code, mask, size, mutfun):
    tsize = random.randint(1,size-1)
    copypos = random.randint(1, len(code)-tsize)
    insertpos = random.randint(1, len(code)-tsize)
    code.insert(code[copypos:copypos+tsize],insertpos)
    mutfun(code,mask)
    return code

def junk(code, mask, size, mutfun):
    tsize = random.randint(1,size-1)
    instream = BitStream(bin='0b'+'0'*tsize)
    insertpos = random.randint(1, len(code)-tsize)
    code.insert(instream,insertpos)
    mutfun(code)
    return code

def delete(code, mask, size, mutfun):
    tsize = random.randint(1,size-1)
    delpos = random.randint(1, len(code)-tsize)
    del code[delpos:delpos+tsize]
    mutfun(code,mask)
    return code

#FIXME: will not work correctly with variable size genomes
def crossover(code1, code2, numpoints = 1,mutrate=0.01):
    offspring1 = BitStream()
    offspring2 = BitStream()
    r = random.random()
    if r < 0.9:
        minsize = min([len(code1),len(code2)])
        points = [random.randint(0,minsize-1) 
                  for p in range(numpoints)]
        points.append(minsize)
        points.sort()
        last = 0
        switch = False;
   
        for p in points:
            if not switch:
                offspring1 += code1[last:p]
                offspring2 += code2[last:p]
            else:
                offspring1 += code2[last:p]
                offspring2 += code1[last:p]
            last = p
            switch = not switch
        else:
            offspring1 = BitStream(code1)
            offspring2 = BitStream(code2)
        

    bitflipmutation(offspring1,mutrate)
    bitflipmutation(offspring2,mutrate)
    return [offspring1, offspring2]
        
    

if __name__ == '__main__':
    a = BitStream('0b' + '1'*20)
    b = BitStream('0b' + '0'*20)
    c = BitStream('0b' + '1'*20)
    print a.bin
    transposon(a, 0.01, 0.5, 20)
    print a.bin
    junk(a, 0.01, 0.5, 20)
    print a.bin
    delete(a, 0.01, 0.5, 20)
    print a.bin
    offsprings = crossover(c,b)
    print offsprings[0].bin   
    print offsprings[1].bin
