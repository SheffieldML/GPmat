import numpy as np

def getSymbols(number):

    """Get a cell array of different plot symbols.
    
    Description:
    
    symbol = getSymbols(number) returns a cell array of different plot
     symbols. A maximum of 66 distinct symbols will be created.
     Returns:
      symbol - cell array of the different symbols.
     Arguments:
      number - the number of different plot symbols required.
        

    See also
    PLOT


    Copyright (c) 2005, 2009 Neil D. Lawrence
    
    """


    symbolColour = ['r', 'g', 'b', 'c', 'm']
    symbolShape = ['x', 'o', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p']
    counter = 0
    symbol = []
    while counter < number:
        symbol.append(symbolColour[counter % len(symbolColour)] + symbolShape[counter % len(symbolShape)])
        counter += 1

    return symbol

def readBinaryDoubles(fileName):
    """READBINARYDOUBLES Read information from a binary file in as doubles.
    FORMAT
    DESC reads in information from a binary file as a vector of doubles.
    ARG fileName : the name of the file.
    RETURN vec : vector of values in the file. 
    
    COPYRIGHT : Neil D. Lawrence, 2009
    
    SEEALSO : open, read, close"""
 

    fid = open(fileName)
    X = np.fromfile(fid, dtype=float, count=-1)
    return X


def  jitChol(A, maxTries=10, warning=True):

    """Do a Cholesky decomposition with jitter.
    
    Description:
    
    
    U, jitter = jitChol(A, maxTries, warning) attempts a Cholesky
     decomposition on the given matrix, if matrix isn't positive
     definite the function adds 'jitter' and tries again. Thereafter
     the amount of jitter is multiplied by 10 each time it is added
     again. This is continued for a maximum of 10 times.  The amount of
     jitter added is returned.
     Returns:
      U - the Cholesky decomposition for the matrix.
      jitter - the amount of jitter that was added to the matrix.
     Arguments:
      A - the matrix for which the Cholesky decomposition is required.
      maxTries - the maximum number of times that jitter is added before
       giving up (default 10).
      warning - whether to give a warning for adding jitter (default is True)

    See also
    CHOL, PDINV, LOGDET


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """

    jitter = 0
    i = 0
    while(True):
        try:
            # Try --- need to check A is positive definite
            if jitter == 0:
                jitter = abs(np.trace(A))/A.shape[0]*1e-6
                LC = np.linalg.cholesky(A)
                return LC.T, 0.0
            else:
                if warning:
                    print "Adding jitter of ", jitter, " in jitChol()."
                LC = np.linalg.cholesky(A+jitter*np.eye(A.shape[0]))
                print 
                return LC.T, jitter
        except np.linalg.LinAlgError:
            # Seems to have been non-positive definite.
            if i<maxTries:
                jitter = jitter*10
            else:
                raise np.linalg.LinAlgError, "Matrix non positive definite, jitter of " +  str(jitter) + " added but failed after " + str(i) + " trials."
        i += 1

def pdinv(A, UC=None, maxTries=10, warning=True):

    """Invert a positive definite matrix.
    
    Description:
    
    [Ainv, U] = pdinv(A) inverts a positive definite matrix. If the
     matrix isn't quite positive definite the function adds 'jitter' to
     make it positive definite and gives out a warning message (this is
     done through JITCHOL).
     Returns:
      Ainv - the inverse of A computed using Cholesky decomposition.
      U - the Cholesky decomposition of A.
     Arguments:
      A - the input positive definite matrix to be inverted.
    
    [Ainv, U] = pdinv(A, U) inverts a positive definite matrix given
     the Cholesky decomposition of A.
     Returns:
      Ainv - the inverse of A computed using Cholesky decomposition.
      U - the Cholesky decomposition of A.
     Arguments:
      A - the input positive definite matrix to be inverted.
      U - the Cholesky decomposition of A.
    
    [Ainv, U, jitter] = pdinv(A, U) inverts a positive definite matrix
     given the Cholesky decomposition of A. If jitter is used then the
     amount of jitter used is returned.
     Returns:
      Ainv - the inverse of A computed using Cholesky decomposition.
      U - the Cholesky decomposition of A.
      jitter - the amount of jitter added.
     Arguments:
      A - the input positive definite matrix to be inverted.
      U - the Cholesky decomposition of A.
        

    See also
    jitChol, logdet, chol


    Copyright (c) 2003, 2004, 2005, 2006, 2009 Neil D. Lawrence
    
    """
    
    # Obtain a Cholesky decomposition.
    if UC==None:
        UC, jitter = jitChol(A, maxTries, warning)
    else:
        jitter = None
    invU = np.asmatrix(np.linalg.solve(UC,np.eye(A.shape[0])))
    
    Ainv = invU*invU.T

    return Ainv, UC, jitter

def logdet(A, UC=None, maxTries=10, warning=True):

    """The log of the determinant when argument is positive definite.
    
    Description:
    
    [d, U] = logdet(A) returns the log determinant of a positive
     definite matrix. If the matrix isn't quite positive definite the
     function adds 'jitter' to make it positive definite and gives out
     a warning message (this is done through JITCHOL).
     Returns:
      d - the log determinant of A computed using Cholesky
       decomposition.
      U - the Cholesky decomposition of A.
     Arguments:
      A - the input positive definite matrix for which the log
       determinant is required.
    
    [d, U, jitter] = logdet(A, U) returns the log determinant of a
     positive definite matrix given the Cholesky decomposition of A. If
     jitter is used then the amount of jitter used is returned.
     Returns:
      d - the log determinant of A computed using Cholesky
       decomposition.
      U - the Cholesky decomposition of A.
      jitter - the amount of jitter added.
     Arguments:
      A - the input positive definite matrix for which the log
       determinant is required.
      U - the Cholesky decomposition of A.
        

    See also
    jitChol, pdinv, chol


    Copyright (c) 2003, 2004, 2005, 2006, 2009 Neil D. Lawrence
    
    """
    
    # Obtain a Cholesky decomposition.
    if UC==None:
        UC, jitter = jitChol(A, maxTries, warning)
        
    ld = 2*np.sum(np.log(np.diag(UC)));
    return ld, UC
