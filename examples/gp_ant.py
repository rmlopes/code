#    This file is part of DEAP.
#
#    DEAP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    DEAP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with DEAP. If not, see <http://www.gnu.org/licenses/>.

"""
This example is from "John R. Koza. Genetic Programming: On the Programming
of Computers by Natural Selection. MIT Press, Cambridge, MA, USA, 1992.".

The problem is called The Artificial Ant Problem.
<http://www.cs.ucl.ac.uk/staff/w.langdon/bloat_csrp-97-29/node2.html>

The goal of this example is to show how to use DEAP and its GP framework with
with complex system of functions and object.

Given an AntSimulator ant, this solution should get the 89 pieces of food
within 543 moves.
ant.routine = ant.if_food_ahead(ant.move_forward, prog3(ant.turn_left,
                                                  prog2(ant.if_food_ahead(ant.move_forward, ant.turn_right),
                                                        prog2(ant.turn_right, prog2(ant.turn_left, ant.turn_right))),
                                                  prog2(ant.if_food_ahead(ant.move_forward, ant.turn_left), ant.move_forward)))

Best solution found with DEAP:
prog3(prog3(move_forward,
            turn_right,
            if_food_ahead(if_food_ahead(prog3(move_forward,
                                              move_forward,
                                              move_forward),
                                        prog2(turn_left,
                                              turn_right)),
                          turn_left)),
      if_food_ahead(turn_left,
                    turn_left),
      if_food_ahead(move_forward,
                    turn_right))
fitness = (89,)
"""

import sys
import logging
import copy
import random
import string
from functools import partial


logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)

def dummy():
    return

def progn(*args):
    #print 'PROGN: ', args
    for arg in args:
        arg()

def prog2(out1, out2):
    return partial(progn,out1,out2)

def prog3(out1, out2, out3):
    return partial(progn,out1,out2,out3)

def progN(*args):
    return partial(progn,*args)

def if_food_ahead_RL(ant, out1, *targs):
        #print targs
        return partial(if_then_else, ant.sense_food, out1, *targs)

#modified for variable arity
#def if_then_else(condition, out1, out2):
    #print targs
 #   if condition(): out1()
  #  else: out2()
     #   for out in targs:
      #      out()

def if_then_else(condition, out1,*targs):
    #print targs
    if condition(): out1()
    else:
       for out in targs:
            out()

class AntSimulator(object):
    direction = ["north","east","south","west"]
    dir_row = [1, 0, -1, 0]
    dir_col = [0, 1, 0, -1]

    def __init__(self, max_moves):
        self.max_moves = max_moves
        self.moves = 0
        self.eaten = 0
        self.routine = None

    def _reset(self):
        self.row = self.row_start
        self.col = self.col_start
        self.dir = 1
        self.moves = 0
        self.eaten = 0
        self.matrix_exc = copy.deepcopy(self.matrix)#self.matrix[:]

    #NOTE: with annotations is not possible to use the functions as
    #      named arguments in progn. With eval the @property is not necessary.
    #For multi-line programs (uses exec) the full notation must be used
    #    eg.: 'if ant.sense_food():\n\tant.move_forward()\n'
    @property
    def position(self):
        return (self.row, self.col, self.direction[self.dir])

    #@property
    def turn_left(self):
        if self.moves < self.max_moves:
            self.moves += 1
            self.dir = (self.dir - 1) % 4
    #@property
    def turn_right(self):
        if self.moves < self.max_moves:
            self.moves += 1
            self.dir = (self.dir + 1) % 4

    #@property
    def move_forward(self):
        if self.moves < self.max_moves:
            self.moves += 1
            self.row = (self.row + self.dir_row[self.dir]) % self.matrix_row
            self.col = (self.col + self.dir_col[self.dir]) % self.matrix_col
            if self.matrix_exc[self.row][self.col] == "food":
                self.eaten += 1
            self.matrix_exc[self.row][self.col] = "passed"

    def sense_food(self):
        ahead_row = (self.row + self.dir_row[self.dir]) % self.matrix_row
        ahead_col = (self.col + self.dir_col[self.dir]) % self.matrix_col
        return self.matrix_exc[ahead_row][ahead_col] == "food"

    #modified for variable arity
    #evaluation with exec requires to inverse here the arguments
    #since they are read from right to left
    #validate in gearnet
    #def if_food_ahead(self, out1, out2):
      #  return partial(if_then_else, self.sense_food, out1, out2)

    def if_food_ahead(self, out1, *targs):
        return partial(if_then_else, self.sense_food, out1, *targs)

    def run(self):
        self._reset()
        while self.moves < self.max_moves:
            try:
                self.routine()
            except SyntaxError:
                print "SYNTAX ERROR:\n"+routine
                exit(0)


            #FIXME: its found now fix it
    def runstring(self,routine,callable_ = False):
        self._reset()
        while self.moves < self.max_moves:
            last = self.moves
            try:
                if callable_:
                    eval(routine,{'ant':self,'prog3':prog3,
                                  'prog2':prog2, 'progN':progN, 'dummy':dummy})()
                else:
                    exec(routine,{'ant':self})

                #protection for circuits that don't make the ant move
                #loosing energy
                if last == self.moves:
                    self.moves+=1
            except SyntaxError:
                print "SYNTAX ERROR:\n"+routine
                exit(0)


    def parse_matrix(self, matrix):
        self.matrix = list()
        for i, line in enumerate(matrix):
            self.matrix.append(list())
            for j, col in enumerate(line):
                if col == "#":
                    self.matrix[-1].append("food")
                elif col == ".":
                    self.matrix[-1].append("empty")
                elif col == "S":
                    self.matrix[-1].append("empty")
                    self.row_start = self.row = i
                    self.col_start = self.col = j
                    self.dir = 1
        self.matrix_row = len(self.matrix)
        self.matrix_col = len(self.matrix[0])
        self.matrix_exc = copy.deepcopy(self.matrix)

ant = AntSimulator(615)
trail_file = open("datafiles/santafe_trail.txt")
ant.parse_matrix(trail_file)

if __name__ == "__main__":
    #routine='progn(progn(progn(ant.turn_left),ant.turn_left),progn(ant.turn_left),ant.turn_left)'
    #routine='progn(ant.move_forward, progn(ant.turn_right, ant.turn_left))'
    #head(ant.move_forward,ant.turn_right,ant.if_food_ahead(ant.move_forward,ant.turn_right,ant.if_food_ahead(ant.move_forward,ant.turn_right,ant.move_forward)))))'
#f = open('ant-success.txt')
    #lines = f.readlines()
    #routine = reduce(lambda l1,l2: l1 + l2, lines)
    #print routine
#    ant.routine = 'ant.if_food_ahead(ant.move_forward, prog3(ant.turn_left, prog2(ant.if_food_ahead(ant.move_forward, ant.turn_right), prog2(ant.turn_right, prog2(ant.turn_left, ant.turn_right))), prog2(ant.if_food_ahead(ant.move_forward, ant.turn_left), ant.move_forward)))'
    #testantroutine(routine)
    #routine = 'ant.if_food_ahead(ant.move_forward, prog3(ant.turn_left,prog2(ant.if_food_ahead(ant.move_forward, ant.turn_right), prog2(ant.turn_right, prog2(ant.turn_left, ant.turn_right))), prog2(ant.if_food_ahead(ant.move_forward, ant.turn_left), ant.move_forward)))'
    #ant.routine = ant.if_food_ahead(
      #  ant.move_forward,
    #    progN(ant.turn_left,
     #         progN(ant.if_food_ahead(ant.move_forward, ant.turn_right),
    #                prog2(ant.turn_right, prog2(ant.turn_left, ant.turn_right))),
      #        progN(ant.if_food_ahead(ant.move_forward, ant.turn_left),
       #             ant.move_forward)))
    #ant.run()

    mystr = 'if ant.sense_food():\n\tant.move_forward()\n'
    mystr2 = 'ant.if_food_ahead(ant.move_forward,ant.move_forward)'
    print mystr
    ant.routine = mystr
    #ant.run()
    ant.runstring(mystr2, True)
    print ant.eaten
