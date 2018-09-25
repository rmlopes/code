from arn import *
import numpy as np

def getmoore(x, y):
    return [(x,y+1), (x+1, y+1), (x+1,y), (x+1, y-1),
            (x,y-1), (x-1, y-1), (x-1,y), (x-1, y+1)]

def getvonneumann(x, y):
        return [(x,y+1), (x+1,y),
                (x,y-1), (x-1,y)]

_CELLSIZE = 5
_INTENSITY = 2
_COEFFICIENT = .2

if __name__ == '__main__':
        arnconfigfile = '../configfiles/arnsim.cfg'
        log.setLevel(logging.DEBUG)
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(open(arnconfigfile))
        proteins=[]
        nump = 0
        try:
                f = open(sys.argv[1], 'r')
                genome = BitStream(bin=f.readline())
                arnet = ARNetwork(genome, cfg)
        except:
                while nump < 3:
                        genome = BitStream(float=random.random(), length=32)
                        for i in range(cfg.getint('default','initdm')):
                                genome = dm_event(genome,
                                                  .02)

                                arnet = ARNetwork(genome, cfg)
                                nump = len(arnet.promlist)

        f = open('genome.save','w')
        f.write(genome.bin)
        f.close

        for p in arnet.proteins: print p
        c = Cell(arnet,_CELLSIZE)

 #   for x in range(_CELLSIZE):
  #      for y in range(_CELLSIZE):
   #         for k in arnet.promlist:
    #            print [c.proteins[k][e[0]][e[1]]
     #                  for e in getneighbors(x,y)+[(x,y)]
      #                 if  0 < e[0] < c.cellsize and
       #                0 < e[1] < c.cellsize]

        cchistory=nparray(c.ccs)
        for i in range(2):
                c.stepsimulate()
                cchistory = np.column_stack((cchistory,c.ccs))

        displayARNresults(arnet.proteins, cchistory)
    #c.diffuseproteins(_COEFFICIENT,_INTENSITY)
    #for k,v in c.proteins.items():
     #   print k, v