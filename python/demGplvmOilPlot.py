import numpy
import matplotlib.pyplot as pyplot
import ndlwrap

def plotOil(X, lbls):
    Xstor = ndlwrap.toarray(X)


    pyplot.figure()
    ind = numpy.nonzero(lbls[:, 0]==1)
    pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'ro')
    ind = numpy.nonzero(lbls[:, 1]==1)
    pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'bx')
    ind = numpy.nonzero(lbls[:, 2]==1)
    pyplot.plot(Xstor[ind, 0], Xstor[ind, 1], 'gs')
