import sys
sys.path.insert(0,'..')
import numpy as np
import pylab as pb
from gp import chgp
from gp import hgp
import kern

pb.ion()

#build a double-hierarchy HGP
#1) construct the data
Nrep = 5
Ngene = 4
Nd = [np.random.randint(2,8) for i in range(Nrep)]
X = [np.random.randn(Ndi,1) for Ndi in Nd]
[xx.sort(0) for xx in X]
Y = [[np.sin(xx) + 0.3*np.sin(xx+5*np.random.rand()) + 0.5*np.sin(xx+10.*gs) + 0.05*np.random.randn(Ndi,1) for Ndi,xx in zip(Nd,X)] for g,gs in enumerate([np.random.randn() for i in range(Ngene)])]

#2) construct the pdata
pdata_chgp = np.vstack([np.tile([[str(i)]],(Ndi,1)) for i,Ndi in enumerate(Nd)])
pdata_rep = np.tile(pdata_chgp,(Ngene,1))
pdata_gene = np.vstack([np.tile([[str(i)]],np.vstack(X).shape) for i in range(Ngene)])
pdata_hgp = np.hstack((pdata_gene,np.array(['_'.join(np.hstack((g,r))) for g,r in zip(pdata_gene, pdata_rep)] )[:,None]))

#3) construct the hgp model
m1 = hgp(np.tile(np.vstack(X),(Ngene,1)),np.vstack([np.vstack(yy) for yy in Y]),pdata_hgp)
m1.constrain_positive('')
m1.optimize()
m1.plot(1)

#4) construct nested chgp model
tmp = hgp(np.vstack(X),np.random.randn(np.vstack(X).shape[0],1),pdata_chgp)
m2 = chgp(np.vstack(X),np.hstack([np.vstack(yy) for yy in Y]),kern.rbf(np.vstack(X)),tmp.kern)
m2.constrain_positive('')
m2.optimize()
m2.plot()


