from math import log,exp,tan,tanh,sqrt,cos,sin,cosh,sinh
import operator

#### MATH ###
def add_(x1, x2 = 0):
    return operator.add(x1,x2)

def sub_(x1, x2 = 0):
    return operator.sub(x1,x2)

def mul_(x1, x2 = 1):
    return operator.mul(x1,x2)

def div_(x, y = 1):
    if(y==0): return 1
    return x / y

def tanh_(*args):
    inp = sum(args)
    return tanh(inp)

def tan_(*args):
    inp = sum(args)
    return tanh(inp)
def sin_(*args):
    inp = sum(args)
    return tanh(inp)
def cos_(*args):
    inp = sum(args)
    return tanh(inp)
def sinh_(*args):
    inp = sum(args)
    return tanh(inp)
def cosh_(*args):
    inp = sum(args)
    return tanh(inp)

def log_(*args):
    x = sum(args)
    if x <= 0: return 0
    return log(x)

def exp_(*args):
    x= sum(args)
    try:
        return exp(x)
    except:
        return 1

def sqrt_(*args):
    x = sum(args)
    return sqrt(x)

### LOGIC ###
def nand(in1, in2 = 1):
    result = not(and_(in1, in2))
    #print 'in1=', in1, ' in2=',in2, ' r=',result
    return result

def nor(in1, in2 = 0):
    result = not(or_(in1, in2))
    #print 'in1=', in1, ' in2=',in2, ' r=',result
    return result

def and_(in1, in2 = 1):
    result = (int(in1) and int(in2)) 
    #print 'in1=', in1, ' in2=',in2, ' r=',result
    return result

def or_(in1, in2 = 0):
    result = (int(in1) or int(in2))
    #print 'in1=', in1, ' in2=',in2, ' r=',result
    return result


### CLASSIFICATION ###

def mcc(tp, tn, fp, fn):
    '''Matthews Correlation Coefficient'''
    den = sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
    if not den:
        den = 1
    return (tp * tn - fp * fn) / den


#FROM http://code.activestate.com/recipes/384122-infix-operators/
# definition of an Infix operator class
# this recipe also works in jython
# calling sequence for the infix is either:
#  x |op| y
# or:
# x <<op>> y
class Infix:
    def __init__(self, function):
        self.function = function

    def __ror__(self, other):
        return Infix(lambda x, self=self, other=other: self.function(other, x))
    def __or__(self, other):
        return self.function(other)
    def __rlshift__(self, other):
        return Infix(lambda x, self=self, other=other: self.function(other, x))
    def __rshift__(self, other):
        return self.function(other)
    def __call__(self, value1, value2):
        return self.function(value1, value2)

#this was used with grammars to evaluate expressions
div = Infix(div_)

if __name__=='__main__':
    print 8 |_div_| 2
    print 9.0 |_div_| 2
    print 8 |_div_| 0
    print 8 / 0
