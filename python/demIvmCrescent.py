#!/usr/bin/env python

# Test the Python/C++ IVM on the crescent data.
import ndlml
import ndlmatrix
import numpy
import ndldata

X, y = ndldata.generateCrescentData(40)

X = ndlmatrix.ndlmatrix(X)
y = ndlmatrix.ndlmatrix(y)

# Set up kernel function
kern = ndlml.cmpndKern()

kern1 = ndlml.rbfKern(X)
kern2 = ndlml.biasKern(X)
kern3 = ndlml.whiteKern(X)

kern.addKern(kern1)
kern.addKern(kern2)
kern.addKern(kern3)

noise = ndlml.probitNoise(y)

# Create an IVM model.
model = ndlml.ivm(X, y, kern, noise, ndlml.ivm.ENTROPY, 100, 3)

# Optimise the IVM.
model.optimise(14, 100, 0)

