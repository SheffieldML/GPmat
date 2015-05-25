#This code is lifted directly from the scipy cookbook, http://www.scipy.org/Cookbook/MultiDot, and remains copyright of scipy.org.

import types
import numpy
def mdot(*args):
   """Multiply all the arguments using matrix product rules.
   The output is equivalent to multiplying the arguments one by one
   from left to right using dot().
   Precedence can be controlled by creating tuples of arguments,
   for instance mdot(a,((b,c),d)) multiplies a (a*((b*c)*d)).
   Note that this means the output of dot(a,b) and mdot(a,b) will differ if
   a or b is a pure tuple of numbers.
   """
   if len(args)==1:
       return args[0]
   elif len(args)==2:
       return _mdot_r(args[0],args[1])
   else:
       return _mdot_r(args[:-1],args[-1])

def _mdot_r(a,b):
   """Recursive helper for mdot"""
   if type(a)==types.TupleType:
       if len(a)>1:
           a = mdot(*a)
       else:
           a = a[0]
   if type(b)==types.TupleType:
       if len(b)>1:
           b = mdot(*b)
       else:
           b = b[0]
   return numpy.dot(a,b)
