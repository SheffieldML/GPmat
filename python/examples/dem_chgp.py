import sys
sys.path.append('..')
import numpy as np
import pylab as pb
from gp import chgp
from gp import hgp
import kern

pb.ion()

#A chgp demo:
Nd = 30
Ng = 8
X = np.random.randn(Nd,1)
X.sort(0)
Y = np.sin(X) + 0.3*np.sin(X+5*np.random.rand(1,Ng)) + np.random.randn(Nd,Ng)*0.1
kernf = kern.rbf(X)
kerny = kern.rbf(X) + kern.white(X)

m1 = chgp(X,Y,kernf,kerny)
m1.constrain_positive('')
m1.optimize()
m1.plot()

#build an equivalent hgp model
XX = np.tile(X,(Ng,1))
YY = Y.T.flatten()[:,None]
pd = np.vstack([np.tile([[str(i)]],(Nd,1)) for i in range(Ng)])
m2 = hgp(XX,YY,pd)
m2.kern.unconstrain('')
m2.constrain_positive('')
m2.optimize()
m2.plot()

