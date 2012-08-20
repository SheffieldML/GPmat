import numpy as np
import pylab as pb
import GPy

def go():
	#generate GPLVM-like data
	X = np.random.randn(50,1)
	k = GPy.kern.rbf(X) + GPy.kern.white(X)
	k.kerns[1].alpha = 0.0001
	k.kerns[0].gamma = 0.002
	K = k.compute()
	Y = np.random.multivariate_normal(np.zeros(50),K,2).T
	
	Y = Y-Y.mean(0)
	Y /= Y.std(0)
	print "y first"
	print Y.shape
	print max(Y[0])
	print min(Y[0])
	print X.shape
	#find PCA solution and construct GPLVM
	XPCA,W = GPy.util.linalg.PCA(Y,2)
	k = GPy.kern.rbf(XPCA) + GPy.kern.white(XPCA)
	k.kerns[1].alpha = np.array([0.001])
	k.kerns[0].gamma = np.array([0.02])
	k.constrain_positive('cmp')
	m = GPy.models.GPLVM(Y,XPCA,2,k)

	#m.checkgrad()

	m.plot()
	pb.title('PCA initialisation')

	m.optimize(max_f_eval=50)
	y = m.plot()
	print min(y[:,0])
	pb.title('After optimisation')


