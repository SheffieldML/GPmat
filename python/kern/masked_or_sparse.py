from timeit import timeit
import sys
import numpy as np

setup =\
"""
import numpy as np
from scipy import sparse
A = np.random.rand(%d,%d)
M = np.random.binomial(1,%f,(%d,%d))
Ama = np.ma.masked_array(A,M)
i1,i2 = np.nonzero(M)
Asp = sparse.coo.coo_matrix((A[i1,i2],(i1,i2)),shape=(%d,%d))

B = np.random.rand(%d,%d)
M = np.random.binomial(1,%f,(%d,%d))
Bma = np.ma.masked_array(B,M)
i1,i2 = np.nonzero(M)
Bsp = sparse.coo.coo_matrix((A[i1,i2],(i1,i2)),shape=(%d,%d))
"""

sp_kern = "np.asarray((Asp+Bsp).todense())"
ma_kern = "np.exp(Ama).filled()+np.exp(Bma).filled()"

res_sp = np.zeros((3,4))
res_ma = np.zeros((3,4))
for n,N in enumerate([10,100,1000]):
	for r,R in enumerate([.01,.1,.2,.5]):
		res_sp[n,r] = timeit(setup=setup%(N,N,R,N,N,N,N, N,N,R,N,N,N,N),stmt=sp_add,number=10000/N)
		res_ma[n,r] = timeit(setup=setup%(N,N,R,N,N,N,N, N,N,R,N,N,N,N),stmt=ma_add,number=10000/N)
		print N, r
		sys.stdout.flush()
