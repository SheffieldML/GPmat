import os
import posix
import numpy
import ndlml
import scipy.sparse
import scipy.io
import matplotlib.mlab

import pdb

baseDir = os.path.join(posix.environ['HOME'], 'mlprojects', 'datasets', 'matlab')


def lvmLoadData(dataSetName, seedVal=10000):

    """This is for loading in data sets for processing using machine
    learning algorithms, particularly latent variable models."""

    numpy.random.seed(seed=seedVal)

    if dataSetName == 'oil':
        fid = open(os.path.join(baseDir, 'oil', 'DataTrn.txt'))
        Y = numpy.fromfile(fid, sep='\t').reshape((-1, 12))
        fid.close()
        fid = open(os.path.join(baseDir, 'oil', 'DataTrnLbls.txt'))
        lbls = numpy.fromfile(fid, sep='\t').reshape((-1, 3))
        fid.close()
        return Y,lbls

    elif dataSetName == 'swissRoll':
        matData = scipy.io.loadmat(os.path.join(baseDir, 'swiss_roll_data'))
        Y = matData['X_data'][:, 0:1000].transpose()
        return Y
    
    elif dataSetName == 'swissRollFull':
        matData = scipy.io.loadmat(os.path.join(baseDir, 'swiss_roll_data'))
        Y = matData['X_data'][:, 0:3000].transpose()
        return Y
    
    elif dataSetName == 'oil100':
        Y, lbls = lvmLoadData('oil')
        indices = numpy.random.permutation(1000)
        indices = indices[0:100]
        Y = Y[indices, :]
        lbls = lbls[indices, :]
        return Y,lbls
    
    elif dataSetName[0:-1] == 'movielensSmall':
        partNo = int(dataSetName[-1])
        fileName = os.path.join(baseDir, 'movielens', 'small', 'u' + str(partNo) + '.base')
        fid = open(fileName)
        uTrain = numpy.fromfile(fid, sep='\t', dtype=numpy.int16).reshape((-1, 4))
        fid.close()
        maxVals = numpy.amax(uTrain, axis=0)
        numUsers = maxVals[0]
        numFilms = maxVals[1]
        numRatings = uTrain.shape[0]
        
        Y = scipy.sparse.lil_matrix((numFilms, numUsers), dtype=numpy.int8)
        for i in range(numUsers):
            ind = matplotlib.mlab.find(uTrain[:, 0]==i+1)
            Y[uTrain[ind, 1]-1, i] = uTrain[ind, 2] 
        
        fileName = os.path.join(baseDir, 'movielens', 'small', 'u' + str(partNo) + '.test')
        fid = open(fileName)
        uTest = numpy.fromfile(fid, sep='\t', dtype=numpy.int16).reshape((-1, 4))
        fid.close()
        numTestRatings = uTest.shape[0]

        Ytest = scipy.sparse.lil_matrix((numFilms, numUsers), dtype=numpy.int8)
        for i in range(numUsers):
            ind = matplotlib.mlab.find(uTest[:, 0]==i+1)
            Ytest[uTest[ind, 1]-1, i] = uTest[ind, 2] 

        lbls = numpy.empty((1,1))
        lblstest = numpy.empty((1,1))
        return Y,lbls,Ytest,lblstest

def mapLoadData(dataSetName, seedVal=100000):

    """This is for loading in data sets for processing using machine
    learning algorithms, particularly mapping models (classification,
    regression) latent variable models.

    Options for dataSetName include 'spgp1d', 'gp1d' """

    numpy.random.seed(seed=seedVal)
    if dataSetName == 'spgp1d':
        numIn = 1
	N = 500
	X = numpy.random.uniform(low=-1.0, high=1.0, size=(N, numIn))
	X.sort(axis=0)
	Xn = ndlml.matrix(X)
	kern1 = ndlml.rbfKern(Xn)
	kern1.setVariance(1)
	kern1.setInverseWidth(20)
	kern2 = ndlml.whiteKern(Xn)
	kern2.setVariance(0.01)
	kern = ndlml.cmpndKern(Xn)
	kern.addKern(kern1)
	kern.addKern(kern2)
	K = ndlml.matrix(numpy.empty((N, N)))
	kern.compute(K, Xn)
	Kmat = K.toarray()
	y = numpy.reshape(numpy.random.multivariate_normal(numpy.zeros(N), Kmat), (N,1))
	return X,y

    elif dataSetName == 'gp1d':
        numIn = 1
	N = 50
	X = numpy.random.uniform(low=-1.0, high=1.0, size=(N, numIn))
	X.sort(axis=0)
	Xn = ndlml.matrix(X)
	kern1 = ndlml.rbfKern(Xn)
	kern1.setVariance(1)
	kern1.setInverseWidth(20)
	kern2 = ndlml.whiteKern(Xn)
	kern2.setVariance(0.0001)
	kern = ndlml.cmpndKern(Xn)
	kern.addKern(kern1)
	kern.addKern(kern2)
	K = ndlml.matrix(numpy.empty((N, N)))
	kern.compute(K, Xn)
	Kmat = K.toarray()
	y = numpy.reshape(numpy.random.multivariate_normal(numpy.zeros(N), Kmat), (N,1))
	return X,y

    else:
        raise NameError, "Unknown data set"
        


def generateCrescentData(numDataPart):
    sqrt2 = numpy.sqrt(2)
    # Generate a toy data-set
    R = numpy.array([[sqrt2/2, -sqrt2/2], [sqrt2/2, sqrt2/2]])
    D = numpy.array([[2, 0], [0, 1]])
    meanOne = numpy.array([4, 4])
    meanTwo = numpy.array([0, 4])
    meanThree = numpy.array([-4, -4])
    meanFour = numpy.array([0, -4])


    X1 = numpy.random.normal(size=(numDataPart, 2))
    X1 = numpy.mat(X1)*numpy.mat(D)*numpy.mat(R) - meanOne

    X2 = numpy.random.normal(size=(numDataPart, 2))
    X2 = numpy.mat(X2)*numpy.mat(D)*numpy.mat(R) - meanTwo

    X3 = numpy.random.normal(size=(numDataPart, 2))
    X3 = numpy.mat(X2)*numpy.mat(D)*numpy.mat(R) - meanThree

    X4 = numpy.random.normal(size=(numDataPart, 2))
    X4 = numpy.mat(X2)*numpy.mat(D)*numpy.mat(R) - meanFour

    
    X = numpy.vstack((X1, X2, X3, X4))
    y = numpy.vstack((numpy.ones((2*numDataPart, 1)), -numpy.ones((2*numDataPart, 1))))
    return X,y
