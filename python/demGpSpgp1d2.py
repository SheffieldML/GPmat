#!/usr/bin/env python

# Test the Python/C++ GP on the 1-D data.

##SETUP
import sys
import os
import posix
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'swig', 'src'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'datasets', 'python'))
##ENDSETUP

import ndlml as nl
import numpy as np
import datasets
import matplotlib.pyplot as plt
import math

X, y = datasets.mapLoadData('spgp1d')

numActive = 9

Xstore = np.array(X)
ystore = np.array(y)

X = nl.matrix(X)
y = nl.matrix(y)

# Set up kernel function
kern = nl.cmpndKern()

kern1 = nl.rbfKern(X)
kern2 = nl.biasKern(X)
kern3 = nl.whiteKern(X)
kern3.setVariance(1e-4)

kern.addKern(kern1)
kern.addKern(kern2)
kern.addKern(kern3)

noise = nl.gaussianNoise(y)

# Create a GP model.
model = nl.gp(1, 1, X, y, kern, noise, nl.gp.FITC, numActive, 3)
model.setBetaVal(math.exp(2))
model.setDefaultOptimiser(nl.gp.CG)

# Optimise the GP.
model.optimise(1000)


xtest = nl.matrix(np.linspace(-1.2, 1.2, 200).reshape(200,1))
ytest = nl.matrix(200, 1)
yvar = nl.matrix(200, 1)
model.posteriorMeanVar(ytest, yvar, xtest)

xt = xtest.toarray()
yt = ytest.toarray()
yv = yvar.toarray()

tuple1 = (xt, xt[ ::-1,:])
tuple2 = (yt-2*np.sqrt(yv), (yt+2*np.sqrt(yv))[ ::-1,:])

plt.figure()
plt.fill(np.vstack(tuple1), np.vstack(tuple2), facecolor='#a8a8a8', edgecolor='none')
plt.plot(xt, yt, 'k-')

plt.plot(Xstore, ystore, 'k.')

X_u = np.empty((numActive, model.getInputDim()))
X_u = model.X_u.toarray()

plt.plot(X_u, -0.9*np.ones((numActive, 1)), 'bx')
plt.show()
