import sys
sys.path.append('..')
import numpy as np
import pylab as pb
from gp import chgp
from gp import hgp
import kern

#build a double-hierarchy HGP
#1) construct the data
Nrep = 5
Ngene = 4
Nd = [np.random.randint(6,14) for i in range(Nrep)]
X = [np.random.randn(Ndi,1) for Ndi in Nd]
[xx.sort(0) for xx in X]
Y = [[np.sin(xx) + 0.3*np.sin(xx+5*np.random.rand()) + 0.5*np.sin(xx+10.*gs) + 0.05*np.random.randn(Ndi,1) for Ndi,xx in zip(Nd,X)] for g,gs in enumerate([np.random.randn() for i in range(Ngene)])]

#2) construct the pdata
pdata_chgp = np.vstack([np.tile([[str(i)]],(Ndi,1)) for i,Ndi in enumerate(Nd)])
pdata_rep = np.tile(pdata_chgp,(Ngene,1))
pdata_gene = np.vstack([np.tile([[str(i)]],np.vstack(X).shape) for i in range(Ngene)])
pdata_hgp = np.hstack((pdata_gene,np.array(['_'.join(np.hstack((g,r))) for g,r in zip(pdata_gene, pdata_rep)] )[:,None]))

#3) construct the hgp model
m3 = hgp(np.tile(np.vstack(X),(Ngene,1)),np.vstack([np.vstack(yy) for yy in Y]),pdata_hgp)
m3.optimize()
m3.plot(1)

#4) construct nested chgp model
tmp = hgp(np.vstack(X),np.random.randn(np.vstack(X).shape[0],1),pdata_chgp)
m4 = chgp(np.vstack(X),np.hstack([np.vstack(yy) for yy in Y]),kern.rbf(np.vstack(X)),tmp.kern)
m4.kernf.constrain_positive('')
m4.optimize()
m4.plot()


