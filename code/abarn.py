from arn import *
import numpy as np

def diffusion(arr, coef, intensity):
    for i in range(intensity):
        shift = round(random.random())
        if not shift:
            shift = -1
        axis = int(round(random.random()))
        arr+=coef*np.roll(arr,
                          shift=shift,
                          axis=axis)
        #TODO:yield arr

    return arr

def normalize2(cellccs, total):
        for p, ccmap in cellccs.items():
                ccmap /= np.sum(ccmap)

def normalize(cellccs, total):
        for p, ccmap in cellccs.items():
                ccmap /= total

class Cell:
    def __init__(self, arn, cellsize):
        self.proteins = dict([(p[0],
                               np.array([0]*cellsize**2,
                                        dtype=float).reshape(cellsize,cellsize))
                              for p in arn.proteins])
        self.ccs = arn.ccs
        initcc = self.ccs[0]
        center = int(round(cellsize/2))
        for k,v in self.proteins.items():
            v[center][center]  = initcc
        self.cellsize = cellsize
        self.arn = arn

    def diffuseproteins(self, coef, intensity):
        for k,v in self.proteins.items():
            self.proteins[k] = diffusion(v, coef, intensity)

    def stepsimulate(self):
        cellsize = self.cellsize
        newproteins = dict([(p[0],
                             np.array([0]*cellsize**2,
                                      dtype=float).reshape(cellsize,cellsize))
                            for p in self.arn.proteins])

        #totalccs = map(lambda t: (t[0],np.sum(t[1])), self.proteins.items())
        #print totalccs

        for i in range(self.cellsize):
            for j in range(self.cellsize):
                cclist = self.getccs(i,j,self.arn.promlist)
                print cclist
                updated = self.arn.stepsimulate(self.arn.proteins,
                                                cclist)
                #print updated
                ccs = dict( zip( self.arn.promlist,updated))
                for k in self.arn.promlist:
                    newproteins[k][i][j] = ccs[k]

        self.proteins = newproteins
        totalccs = map(lambda t: (t[0],np.sum(t[1])), self.proteins.items())
        #print totalccs
        #
        normalize(self.proteins, sum(zip(*totalccs)[1]))
        totalccs = map(lambda t: (t[0],np.sum(t[1])), self.proteins.items())
        totalccs.sort(key=lambda t: t[0])
        self.ccs = zip(*totalccs)[1]
        #exit(0)
        #for k,v in self.proteins.items():
         #       print k
          #      print v
        print self.ccs

    def getccs(self,x,y, plist):
        result = []
        for k in plist:
            result.append(sum([self.proteins[k][e[0]][e[1]]
                               for e in getneighbors2(x,y)+[(x,y)]
                               if  0 < e[0] < self.cellsize and
                               0 < e[1] < self.cellsize]))
        return result

def getneighbors(x, y):
    return [(x,y+1), (x+1, y+1), (x+1,y), (x+1, y-1),
            (x,y-1), (x-1, y-1), (x-1,y), (x-1, y+1)]

def getneighbors2(x, y):
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