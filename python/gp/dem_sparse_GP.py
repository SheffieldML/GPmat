import numpy as np
import pylab as pb
import sys
import kern
from ndlutil import model
from ndlutil.utilities import pdinv, gpplot
from scipy import optimize
from sparse_GP import vsGP
from simple_GP import GP

#run /home/ricardo/mlprojects/gp/python/gp/dem_sparse_GP.py 
N = 100
M = 10
X = np.random.randn(N,1)
X.sort(0)
Y = np.sin(X)+.5*np.cos(4*X)+np.random.randn(N,1)*0.25
m = vsGP(X,Y,M)
m.constrain_positive('alp|tau|gam')
m.tau = 1.
m.expand_param(m.extract_param())
def f(X):
	m.expand_param(X)
	return  -m.log_likelihood()

m2 = GP(X,Y)
m2.constrain_positive('')

pb.ion()
pb.close('all')
#m.plot()
#m2.plot()

pb.figure()
Xnew = np.linspace(X.min(),X.max(),100)[:,None]
mean_sparse,var_sparse = m.predict(Xnew)
mean_simple,var_simple = m2.predict(Xnew)
gpplot(Xnew,mean_sparse,var_sparse,edgecol='#6495ed',fillcol='#6495ed',alpha=0.6,label='95% CI sparse') 
gpplot(Xnew,mean_simple,var_simple,edgecol='#c0c0c0',fillcol='#c0c0c0',alpha=0.6,label='95% CI full') 
pb.plot(Xnew,mean_sparse,color='red',linewidth=2,label = 'Sparse model')
pb.plot(Xnew,mean_simple,'k--',linewidth=2,label = 'Full model')
pb.plot(X,Y,'kx',mew=1)
pb.plot(m.inducing_inputs(),m.inducing_inputs()*0+pb.ylim()[0],'k|',mew=1.5,markersize=12,label = 'Inducing inputs')
pb.legend(loc='upper left',ncol=2)








