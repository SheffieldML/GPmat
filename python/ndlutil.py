def H = centeringMatrix(dim):

    """returns the centering matrix for the given dimensionality.
    
    Description:
    
    H = centeringMatrix(dim) returns the centering matrix for the
     given dimensionality.
     Returns:
      H - the centering matrix for the given dimensionality.
     Arguments:
      dim - the dimensionality of the centering matrix.
        

    See also
    eye


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    H = -repmat(1./dim, dim, dim) + speye(dim);
def y = cumGamma(x, a, b):

    """Cumulative distribution for gamma.
    
    Description:
    
    p = cumGamma(x) computes the cumulative gamma distribution.
     Returns:
      p - output probability.
     Arguments:
      x - input value.
        

    See also
    gammainc, gamma


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    y = gammainc(x*b, a);

def y = cumGaussian(x):

    """Cumulative distribution for Gaussian.
    
    Description:
    
    p = cumGaussian(x) computes the cumulative Gaussian distribution.
     Returns:
      p - output probability.
     Arguments:
      x - input value.
        

    See also
    lnCumGaussian, lnDiffCumGaussian, erf


    Copyright (c) 2004 Neil D. Lawrence
    
    """
        
    y = 0.5*(1+erf(sqrt(2)/2*x));

def options = defaultOptions;

    """The default options for optimisation.
    
    Description:
    
    options = defaultOptions returns a default options vector for
     optimisation.
     Returns:
      options - the default options vector.
        

    See also
    scg, conjgrad, quasinew


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    options = [0,  1e-4, 1e-4, 1e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-8, 0.1, 0];
def theta = deg2rad(omega):

    """Transform degrees to radians.
    
    Description:
    
    deg2rad(omega, theta) onverts degrees to radians.
     Arguments:
      omega - angle in degrees.
      theta - angle in radians.


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    theta = omega/180*pi;
def [p, lnp] = gammaPdf(x, a, b):

    """PDF for the Gamma distribution.
    
    Description:
    
    [p, lnp] = gammaPdf(x, a, b) computes the pdf of the gamma
     distribution.
     Returns:
      p - probability of the gamma distribution.
      lnp - log of the probability.
     Arguments:
      x - locations where the PDF is to be computed.
      a - shape parameter of the gamma distribution.
      b - rate parameter of the gamma distribuion (inverse scale).
        

    See also
    % SEEALSO cumGamma


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    
    lnp = a*log(b) - gammaln(a) + (a-1)*log(x) - b*x;
    p = exp(lnp);
def y = gaussOverDiffCumGaussian(x, xp, order):

    """A Gaussian over difference of cumulative Gaussians.
    
    Description:
    
    y = gaussOverDiffCumGaussian(X1, X2, order) computes a Gaussian in
     x divided by the difference between two cumulative Gaussian
     distributions.
     Returns:
      y - returns y = ngaussian(X1)/(cumGaussian(X1)-cumGaussian(X2)) if
       order == 1 and ngaussian(X2)/(cumGaussian(X1)-cumGaussian(X2)) if
       order == 2
     Arguments:
      X1 - the argument of the first, positive, cumulative Gaussian.
      X2 - the argument of the second, negative, cumulative Gaussian.
      order - set to 1 or 2, determines whether X1 or X2 is used in the
       argument of the Gaussian term.
        Calculating this function naively causes problems at extreme values.
        
        

    See also
    lnCumGaussian, erfcx, lnDiffCumGaussian, cumGaussian


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    
    %.5*erfcx(-sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(x) ...
    robustAdd = 1e-300;
    fact = sqrt(2)/2;
    xp2 = xp.*xp;
    x2 = x.*x;
    y = zeros(size(xp));
    switch order
     case 1
      expRatio = exp(.5*(x2-xp2));
      index = find(x<=0);
      y(index) = 2./(erfcx(-fact*x(index)) ...
                     -expRatio(index).*erfcx(-fact*xp(index))+robustAdd); 
      x(index) = NaN;
      
      %.5*erfcx(sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(-x)=exp(.5*x*x)*(1-cumGaussian(x))
      index = find(x>0);
      y(index) = 2./(expRatio(index).*erfcx(fact*xp(index))-erfcx(fact*x(index))+robustAdd); 
      y = y/sqrt(2*pi);
     case 2
      expRatio = exp(.5*(xp2-x2));
      index = find(x<=0);
      y(index) = 2./(expRatio(index).*erfcx(-fact*x(index)) ...
                     -erfcx(-fact*xp(index))+robustAdd); 
      x(index) = NaN;
      
      %.5*erfcx(sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(-x)=exp(.5*x*x)*(1-cumGaussian(x))
      index = find(x>0);
      y(index) = 2./(erfcx(fact*xp(index))-expRatio(index).*erfcx(fact*x(index))+robustAdd); 
      y = y*1/sqrt(2*pi);
     otherwise
      error('Incorrect order, should be set to 1 or 2.')
    end
def y = gaussSamp(Sigma, numSamps):

    """Sample from a Gaussian with a given covariance.
    
    Description:
    
    y = gaussSamp(Sigma, numSamps) samples a given number of samples
     from a Gaussian with a given covariance matrix.
     Returns:
      y - the samples from the Gaussian
     Arguments:
      Sigma - the covariance of the Gaussian to sample from.
      numSamps - the number of samples to take from Gaussian.
        

    See also
    randn, eig


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    [U, V] = eig(Sigma);
    dims = size(Sigma, 1);
    y = randn(numSamps, dims);
    y = y*diag(sqrt(diag(V)));
    y = y*U';
def symbol = getSymbols(number):

    """Get a cell array of different plot symbols.
    
    Description:
    
    symbol = getSymbols(number) returns a cell array of different plot
     symbols. A maximum of 66 distinct symbols will be created.
     Returns:
      symbol - cell array of the different symbols.
     Arguments:
      number - the number of different plot symbols required.
        

    See also
    plot


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    symbolColour = {'r', 'g', 'b', 'c', 'm'}; %, 'y'};
    symbolShape = {'x', 'o', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p'};
    counter = 0;
    while counter < number
      symbol{counter+1} = [symbolColour{rem(counter, length(symbolColour))+1} ...
                          symbolShape{rem(counter, length(symbolShape))+1}];
      counter = counter +1;
    end
def lineStr = getline(FID, comment):

    """Get a line from a file.
    
    Description:
    
    lineStr = getline(fid, comment) gets a line from a file, but
     ignore it if it starts with a comment character.
     Returns:
      lineStr - the next line in the file.
     Arguments:
      fid - the identity of the file to load from.
      comment - the character that indicates a line is a comment
       (default #).
        

    See also
    fgetl, fopen


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    
    if nargin < 2
      comment = '#';
    end
    if length(comment)~=1
      error('Comment can only be one character in length');
    end
    lineStr = fgetl(FID);
    if length(lineStr)==0
      return
    end
    while lineStr(1)==comment
      lineStr = fgetl(FID);
    end
def [dlnPart, m] = gradLnDiffErfs(x1, x2, fact1, fact2):,

    """Compute the gradient of the log difference of two erfs.
    
    Description:
        
        


    Copyright (c) 2007 Antti Honkela
    
    """
        
    m = min(x1.^2, x2.^2);
    dlnPart = 2/sqrt(pi) * (exp(-x1.^2 + m) .* fact1 ...
    			- exp(-x2.^2 + m) .* fact2);

def y = gradLogCumGaussian(x):

    """Gradient of the log of the cumulative Gaussian.
    
    Description:
    
    gradLogCumGaussian(x, y) returns the gradient of the logarithm of
     the cumulative Gaussian. Theoretically this is simply
     ngaussian(x)/cumGaussian(x) but there are problems in the tails of
     the distribution. This function attempts to deal with these
     problems without numerical error creeping in.
     Arguments:
      x - the input to the function.
      y - the gradient of the log cumulative Gaussian.
        

    See also
    ngaussian, cumGaussian, erfcx


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    
    y = zeros(size(x));
    index = find(x>0);
    y(index) = ngaussian(x(index))./cumGaussian(x(index));
    index = find(x<=0);
    y(index) = 1./(sqrt(2*pi)*0.5*erfcx(-sqrt(2)/2*x(index)));

def gradientCheck(params, objectiveFunction, gradientFunction, varargin):

    """Check gradients of objective function.
    
    Description:
    
    gradientCheck(params, objectiveFunction, gradientFunction, ...)
     checks the supplied gradient function and the supplied objective
     function to ensure that the numerical gradients (as computed with
     the objective function) match the analytically computed gradients.
     Arguments:
      params - the parameters at which the gradients will be checked.
      objectiveFunction - function handle for the objective function.
      gradientFunction - function handle for the objective function.
      ... - additional arguments that are passed to the objective and
       gradient functions (after the parameter vector which is always
       assumed to be the first argument passed).
        

    See also
    modelObjective, modelGradient, modelCreate


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    
    if isstr(objectiveFunction)
      objectiveFunction = str2func(objectiveFunction);
    end
    if isstr(gradientFunction)
      gradientFunction = str2func(gradientFunction);
    end
    L = 0;
    change = 1e-6;
    origParams = params;
    for i = 1:length(params)
      params(i) = origParams(i) + change;
      Lplus = objectiveFunction(params, varargin{:});
      params(i) = origParams(i) - change;
      Lminus = objectiveFunction(params, varargin{:});
      diff(i) = (Lplus - Lminus)/(2*change);
      params(i) = origParams(i);
    end
    
    anal = gradientFunction(origParams, varargin{:});
    
    delta = anal -diff;
    if max(abs(delta./max([abs(anal) ones(size(anal))], [], 2)))<1e-4
      fprintf('Gradient Check Passed.\n');
    else
      fprintf('Gradient Check Failed.\n');
      fprintf('       \t\t\tAnaly\t\t\tDiffs\t\t\tDelta\n', i, anal(i), diff(i), delta(i));
    
      for i = 1:length(params)
        if(abs(delta(i)/max([abs(anal(i)) 1]))>=1e-4)
          fprintf('Param %d:\t\t%6.4e\t\t%6.4e\t\t%6.4e\n', i, anal(i), diff(i), delta(i));
        end
      end
    end
def hessianCheck(params, objectiveFunction, hessianFunction, varargin):

    """Check Hessian of objective function.
    
    Description:
    
    hessianCheck(params, objectiveFunction, hessianFunction, ...)
     checks the supplied Hessian function and the supplied objective
     function to ensure that the numerical Hessian (as computed with
     the objective function) match the analytically computed Hessian.
     Arguments:
      params - the parameters at which the Hessian will be checked.
      objectiveFunction - function handle for the objective function.
      hessianFunction - function handle for the objective function.
      ... - additional arguments that are passed to the objective and
       Hessian functions (after the parameter vector which is always
       assumed to be the first argument passed).
        

    See also
    modelObjective, modelHessian, modelCreate


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    if isstr(objectiveFunction)
      objectiveFunction = str2func(objectiveFunction);
    end
    if isstr(hessianFunction)
      hessianFunction = str2func(hessianFunction);
    end
    L = 0;
    change = 1e-4;
    origParams = params;
    for i = 1:length(params)
      for j = 1:length(params)
        if i == j
          params(i) = origParams(i) - 2*change;
          LminusMinus = objectiveFunction(params, varargin{:});
          params(i) = origParams(i) + 2*change;
          LplusPlus = objectiveFunction(params, varargin{:});
          params(i) = origParams(i);
          Lother = objectiveFunction(params, varargin{:});
          diffH(i, j) = (LminusMinus -2*Lother + LplusPlus)...
                  /(4*change*change);
        else
          params(i) = origParams(i) - change;
          params(j) = origParams(j) - change;
          LminusMinus = objectiveFunction(params, varargin{:});
          params(i) = origParams(i) - change;
          params(j) = origParams(j) + change;
          LminusPlus = objectiveFunction(params, varargin{:});
          params(i) = origParams(i) + change;
          params(j) = origParams(j) - change;    
          LplusMinus = objectiveFunction(params, varargin{:});
          params(i) = origParams(i) + change;
          params(j) = origParams(j) + change;    
          LplusPlus = objectiveFunction(params, varargin{:});
          diffH(i, j) = (LminusMinus - LminusPlus - LplusMinus + LplusPlus)...
                  /(4*change*change);
        end
        params(i) = origParams(i);
        params(j) = origParams(j);
      end
    end
    
    H = hessianFunction(origParams, varargin{:});
    
    delta = H -diffH;
    if max(abs(delta./max(max([abs(H) ones(size(H))], [], 2))))<1e-4
      fprintf('Hessian Check Passed.\n');
    else
      fprintf('Hessian Check Failed.\n');
      fprintf('       \t\t\tAnaly\t\t\tDiffs\t\t\tDelta\n');
      
      for i = 1:length(params)
        for j = 1:length(params)
          if(abs(delta(i, j)/max([abs(H(i, j)) 1]))>=1e-4)
            fprintf('Param %d, %d:\t\t%6.4e\t\t%6.4e\t\t%6.4e\n', i, j, H(i,j), ...
                    diffH(i,j), delta(i, j));
          end
        end
      end
    end
def y = invCumGaussian(x):

    """Computes inverse of the cumulative Gaussian.
    
    Description:
    
    y = invCumGaussian(x) computes the inverse of the cumulative
     Gaussian.
     Returns:
      y - the inverse of the cumulative Gaussian.
     Arguments:
      x - value between 0 and 1 to map onto the real line.
        

    See also
    cumGaussian, erfinv


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    y = erfinv(x*2 - 1)*2/sqrt(2);
def y=invSigmoid(x):

    """The inverse of the sigmoid function.
    
    Description:
    
    invSigmoid(x, y) returns the inverse of the sigmoid function
     (which takes the form y=log(x/(1-x)).
     Arguments:
      x - the input to the inverse of the sigmoid (should be between 0
       and 1).
      y - the inverse of the sigmoid.
        

    See also
    sigmoid


    Copyright (c) 2004 Neil D. Lawrence
    
    """
        
    y = log(x./(1-x));
def answer = isoctave

    """Returns true if the software running is Octave.
    
    Description:
    
    answer = isoctave tests if the software is octave or not.
     Returns:
      answer - true if the software is octave.


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    try 
      v = OCTAVE_VERSION;
      answer = true;
    catch
      answer = false;
      return
    end
def [UC, jitter] = jitChol(A, maxTries):

    """Do a Cholesky decomposition with jitter.
    
    Description:
    
    U = jitChol(A, maxTries) attempts a Cholesky decomposition on the
     given matrix, if matrix isn't positive definite the function gives
     a warning, adds 'jitter' and tries again. At the first attempt the
     amount of jitter added is 1e-6 times the mean of the diagonal.
     Thereafter the amount of jitter is multiplied by 10 each time it
     is added again. This is continued for a maximum of 10 times.
     Returns:
      U - the Cholesky decomposition for the matrix.
     Arguments:
      A - the matrix for which the Cholesky decomposition is required.
      maxTries - the maximum number of times that jitter is added before
       giving up (default 10).
    
    [U, jitter] = jitChol(A, maxTries) attempts a Cholesky
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
        

    See also
    chol, pdinv, logdet


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    if nargin < 2
      maxTries = 10;
    end
    jitter = 0;
    for i = 1:maxTries
      try
        % Try --- need to check A is positive definite
        if jitter == 0;
          jitter = abs(mean(diag(A)))*1e-6;
          UC = chol(A);
          break
        else
          if nargout < 2
            warning(['Matrix is not positive definite in jitChol, adding ' num2str(jitter) ' jitter.'])
          end
          UC = chol(real(A+jitter*eye(size(A, 1))));
          break
        end
      catch
        % Was the error due to not positive definite?
        nonPosDef = 0;
        verString = version;
        if str2double(verString(1:3)) > 6.1
          [void, errid] = lasterr;
          if strcmp(errid, 'MATLAB:posdef')
            nonPosDef = 1;
          end
        else
          errMsg = lasterr;
          if findstr(errMsg, 'positive definite')
            nonPosDef = 1;
          end
        end
      end
      if nonPosDef
        jitter = jitter*10;
        if i==maxTries
          error(['Matrix is non positive definite tried ' num2str(i) ...
                 ' times adding jitter, but failed with jitter ' ...
                 'of ' num2str(jitter) '. Increase max tries'])
        end
      else
        error(lasterr)
      end
    end
    
    

def kl = kldivGaussian(mean1, cov1, mean2, cov2):

    """Give the KL divergence between two Gaussians.
    
    Description:
    
    kldivGaussian(mean1, cov1, mean2, cov2) returns the
     Kullback-Leibler divergence between two Gaussians with given means
     and covariances.
     Arguments:
      mean1 - mean of the first Gaussian.
      cov1 - covariance of the first Gaussian.
      mean2 - mean of the second Gaussian.
      cov2 - covariance of the second Gaussian.
        

    See also
    logdet, pdinv


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    [invCov2, U] = pdinv(cov2);
    logDet2 = logdet(cov2, U);
    logDet1 = logdet(cov1);
    N = size(cov1, 1);
    meanDiff = mean1 - mean2;
    if size(meanDiff, 1) == 1
      meanDiff = meanDiff';
    end
    kl = -0.5*(logDet1 - logDet2 - trace(cov1*invCov2) ...
               + N - meanDiff'*invCov2*meanDiff);
    

def y = lnCumGaussSum(u1, u2, w1, w2):

    """The log of the weighted sum of two cumulative Gaussians.
    
    Description:
    
    lnCumGaussSum(u1, u2, w1, w2) returns the logarithm of the
     weighted sum of two cumulative Gaussians.
     Arguments:
      u1 - argument of the first cumulative Gaussian.
      u2 - argument of the second cumulative Gaussian.
      w1 - weight of the first cumulative Gaussian.
      w2 - weight of the second cumulative Gaussian.
        

    See also
    cumGaussian, lnCumGaussian, lnDiffCumGaussian


    Copyright (c) 2004 Neil D. Lawrence
    
    """
        
    y = zeros(size(u1));
    safeCond = u1 > 0 & u2 > 0;
    index = find(safeCond);
    if ~isempty(index)
      y(index) = log(w1.*cumGaussian(u1(index)) ...
                     + w2*cumGaussian(u2(index)));
    end
    index = find(~safeCond & u1>u2);
    if ~isempty(index)
      y(index) = log(w1) + lnCumGaussian(u1(index))...
          + log(1 + w2/w1*exp(lnCumGaussian(u2(index))...
                              -lnCumGaussian(u1(index))));
    end
    index = find(~safeCond & u2>=u1);
    if ~isempty(index)
      y(index) = log(w2) + lnCumGaussian(u2(index))...
          + log(1 + w1/w2*exp(lnCumGaussian(u1(index))...
                              -lnCumGaussian(u2(index))));
    end
def y = lnCumGaussian(x):

    """log cumulative distribution for the normalised Gaussian.
    
    Description:
    
    y = lnCumGaussian(X) computes the logarithm of the cumulative
     Gaussian distribution.
     Returns:
      y - log probability of the value under the cumulative Gaussian.
     Arguments:
      X - input position.
        

    See also
    erf, erfcx, cumGaussian, lnDiffCumGaussian, gaussOverDiffCumGaussian


    Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
    
    """
        
    index = find(x< 0);
    if length(index)
      y(index) = -.5*x(index).*x(index) + log(.5) + log(erfcx(-sqrt(2)/2* ...
                                                        x(index)));
    end
    index = find(x>=0);
    if length(index)
      y(index) = log(cumGaussian(x(index)));
    end
    y=reshape(y, size(x));
def f = lnDiffCumGaussian(u, uprime):

    """Log of the difference between two cumulative Gaussians.
    
    Description:
    
    f = lnDiffCumGaussian(u1, u2) computes the logarithm of the
     difference between two cumulative Gaussian distributions.
     Returns:
      f - where f = log(cumGaussian(u1) - cumGaussian(u2)).
     Arguments:
      u1 - the argument of the first (positive) cumulative Gaussian.
      u2 - the argument of the second (negative) cumulative Gaussian.
        

    See also
    cumGaussian, gaussOverDiffCumGaussian, lnCumGaussian


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    f = log(gaussOverDiffCumGaussian(u, uprime, 1)+1e-300) ...
        + .5*u.*u + .5*log(2*pi);
    f=-f;
def [v, signs] = lnDiffErfs(x1, x2):,

    """Helper function for computing the log of difference
    
    Description:
        of two erfs.
    
    [v, s] = lnDiffErfs(x1, x2) computes the log of the difference of
     two erfs in a numerically stable manner.
     Returns:
      v - log(abs(erf(x1) - erf(x2)))
      s - sign(erf(x1) - erf(x2))
     Arguments:
      x1 - argument of the positive erf
      x2 - argument of the negative erf
    
    v = lnDiffErfs(x1, x2) computes the log of the difference of two
     erfs in a numerically stable manner.
     Returns:
      v - log(erf(x1) - erf(x2))     (Can be complex)
     Arguments:
      x1 - argument of the positive erf
      x2 - argument of the negative erf
        

    See also
    gradLnDiffErfs


    Copyright (c) 2007, 2008 Antti Honkela
    
    """
        
    x1 = real(x1);
    x2 = real(x2);
    
    v = zeros(max(size(x1), size(x2)));
    
    if numel(x1) == 1,
      x1 = x1 * ones(size(x2));
    end
    
    if numel(x2) == 1,
      x2 = x2 * ones(size(x1));
    end
    
    signs = sign(x1 - x2);
    I = signs == -1;
    swap = x1(I);
    x1(I) = x2(I);
    x2(I) = swap;
    
    % Case 1: arguments of different signs, no problems with loss of accuracy
    I1 = sign(x1) ~= sign(x2);
    % Case 2: both arguments are positive
    I2 = (x1 > 0) & (x1 > x2) & ~I1;
    % Case 3: both arguments are negative
    I3 = ~I1 & ~I2;
    
    warnState = warning('query', 'MATLAB:log:logOfZero');
    warning('off', 'MATLAB:log:logOfZero');
    v(I1) = log( erf(x1(I1)) - erf(x2(I1)) );
    v(I2) = log(erfcx(  x2(I2)) ...
    	    - erfcx(x1(I2)) .* exp(x2(I2).^2 - x1(I2).^2)) ...
    	- x2(I2).^2;
    v(I3) = log(erfcx(  -x1(I3)) ...
    	    - erfcx(-x2(I3)) .* exp(x1(I3).^2 - x2(I3).^2)) ...
    	- x1(I3).^2;
    warning(warnState.state, 'MATLAB:log:logOfZero');
    
    if nargout < 2,
      v(I) = v(I) + pi*1i;
    end

def [ld, UC] = logdet(A, UC):

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


    Copyright (c) 2003, 2004, 2005, 2006 Neil D. Lawrence
    
    """
        
    if nargin < 2
      UC=[];
    end
    
    % Obtain a Cholesky decomposition.
    if isempty(UC)
      if nargout > 2
        [UC, jitter] = jitChol(A);
      else
        UC = jitChol(A);
      end
    end
    
    ld = 2*sum(log(diag(UC)));


def y = negLogLogit(x):

    """Function which returns the negative log of the logistic function.
    
    Description:
    
    y = negLogLogit(x) computes the negative log of the logistic
     (sigmoid) function, which is also the integral of the sigmoid
     function.
     Returns:
      y - the negative log of the logistic.
     Arguments:
      x - input locations.
        

    See also
    sigmoid


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    y = log(1+exp(x));
def y = ngaussian(x):

    """Compute a Gaussian with mean 0 and variance 1.
    
    Description:
    
    y = ngaussian(x) computes a the likelihood of a normalised
     Gaussian distribution, i.e. with mean 0 and variance 1.
     Returns:
      y - probability of the input values under the Gaussian.
     Arguments:
      x - input value(s) for which to compute the distribution.
        

    See also
    cumGaussian, gaussOverDiffcumGaussian


    Copyright (c) 2004 Neil D. Lawrence
    
    """
        
    x2 = x.*x;
    y = exp(-.5*x2);
    y = y/sqrt(2*pi);

def str = numsf2str(num, sf):

    """Convert number to a string with a number of significant digits.
    
    Description:
    
    str = numsf2str(num, sf) converts a number to a string with a
     given number of significant digits.
     Returns:
      str - the string containing the number with the given number of
       significant digits.
     Arguments:
      num - the number to convert.
      sf - the number of significant figures to show in the string.
        

    See also
    num2str, fprintf


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    val = chop(num, sf);
    str = num2str(val, sf);
    tail = [];
    ePos = find(str == 'e');
    if ~isempty(ePos)
      tail = str(ePos:end);
      str = str(1:ePos-1);
    end
    ind = 1;
    while str(ind) == '0' | str(ind) == '.'
      ind = ind+1;
    end
    count = 0;
    while(ind <= length(str))
      if str(ind) ~= '.'
        count = count + 1;
      end
      ind = ind +1;
    end
    while count < sf
      str = [str '0'];
      count = count + 1;
    end
    str = [str tail];
def [Ainv, UC, jitter] = pdinv(A, UC):

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


    Copyright (c) 2003, 2004, 2005, 2006 Neil D. Lawrence
    
    """
        
    
    if nargin < 2
      UC=[];
    end
    
    % Obtain a Cholesky decomposition.
    if isempty(UC)
      if nargout > 2
        [UC, jitter] = jitChol(A);
      else
        UC = jitChol(A);
      end
    end
    
    invU = UC\eye(size(A, 1));
    %invU = eye(size(A, 1))/UC;
    Ainv = invU*invU'; 

def preparePlot(limitVal, ax):

    """Helper function for tidying up the plot before printing.
    
    Description:
    
    preparePlot(limitVal, ax) is a helper function for tidying up a
     plot before printing.
     Arguments:
      limitVal - the limits to be applied to the axes.
      ax - the axes to apply the plot preparation to.
        

    See also
    zeroAxes


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    axis equal
    
    if nargin < 2
      ax = gca;
    end
    set(ax, 'xlim', limitVal);
    set(ax, 'ylim', limitVal);
    xlim = get(ax, 'xlim');
    ylim = get(ax, 'ylim');
    
    xlim(1) = floor(xlim(1));
    xlim(2) = ceil(xlim(2));
    ylim(1) = floor(ylim(1));
    ylim(2) = ceil(ylim(2));
    
    set(ax, 'xlim', xlim);
    set(ax, 'ylim', ylim);
    
    dirs(1) = max([xlim(1) ylim(1)]);
    dirs(2) = min([xlim(2) ylim(2)]);
    
    line(dirs, dirs);
    line(dirs, dirs+1);
    line(dirs, dirs-1);
    
    zeroaxes(ax, 0.02, 18, 'times')
    

def printPlot(fileName, directory, directoryHtml): 

    """Print a plot to eps and png files.
    
    Description:
    
    printPlot(fileName, directory, directoryPng) prints a plot to the
     specified file name and directory.
     Arguments:
      fileName - the base name of the file to use.
      directory - the directory to place the eps files in.
      directoryPng - the directory to place png the file in.
        

    See also
    preparePlot


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    if nargin < 2 
      directory = '.';
    end
    if nargin < 3
      png = false;
    else
      png = true;
    end
    fprintf('Printing eps plot ...\n');
    print('-depsc', [directory filesep fileName])
    cmap = colormap;
    if png
      fprintf('Printing png plot ...\n');
      % make smaller for PNG plot.
      pos = get(gcf, 'paperposition');
      origpos = pos;
      pos(3) = pos(3)/2;
      pos(4) = pos(4)/2;
      set(gcf, 'paperposition', pos);
      fontsize = get(gca, 'fontsize');
      set(gca, 'fontsize', fontsize/2);
      lineWidth = get(gca, 'lineWidth');
      set(gca, 'lineWidth', lineWidth*2);
      print('-dpng', [directoryHtml filesep fileName]);
      set(gcf, 'paperposition', origpos);
      set(gca, 'fontsize', fontsize);
      set(gca, 'lineWidth', lineWidth);
    end
    
    fprintf('Printing black and white eps plot ...\n');
    colormap gray
    print('-deps', [directory filesep fileName 'NoColour'])
    colormap(cmap);
def vec = readBinaryDoubles(fileName, format):

    """Read information from a binary file in as doubles.
    
    Description:
    
    vec = readBinaryDoubles(fileName, format) reads in information
     from a binary file as a vector of doubles.
     Returns:
      vec - vector of values in the file.
     Arguments:
      fileName - the name of the file.
      format - the file format for 'fopen', defaults to 'ieee-le'.
        

    See also
    fopen, fread, fclose


    Copyright (c) 2009 Neil D. Lawrence
    
    """
          
      if nargin < 2
        format = 'ieee-le';
      end
      fid = fopen(fileName, 'r', format);
      vec = fread(fid, inf, 'double')';
      fclose(fid);
    end
def val = readBoolFromFID(FID, string):

    """Read a boolean from an FID.
    
    Description:
    
    val = readBoolFromFID(FID, name) reads a boolean from a stream.
     Returns:
      val - value of variable in file.
     Arguments:
      FID - stream to read from.
      name - name of boolean.
        

    See also
    writeBoolToFID, readIntFromFID, readDoubleFromFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
      
    val = str2num(readStringFromFID(FID, string));

def val = readDoubleFromFID(FID, string):

    """Read a double from an FID.
    
    Description:
    
    val = readDoubleFromFID(FID, name) reads a double from a stream.
     Returns:
      val - value of variable in file.
     Arguments:
      FID - stream to read from.
      name - name of double.
        

    See also
    writeDoubleToFID, readIntFromFID, readBoolFromFID, readStringFromFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    val = str2num(readStringFromFID(FID, string));

def val = readIntFromFID(FID, string):

    """Read an integer from an FID.
    
    Description:
    
    val = readIntFromFID(FID, name) reads an integer from a stream.
     Returns:
      val - value of variable in file.
     Arguments:
      FID - stream to read from.
      name - name of integer.
        readBoolFromFID, readStringFromFID
        

    See also
    writeIntToFID, readDoubleFromFID, readIntFromFID, 


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
      
    val = str2num(readStringFromFID(FID, string));

def val = readStringFromFID(FID, string):

    """Read an boolean from an FID.
    
    Description:
    
    val = readStringFromFID(FID, name) reads a string from a stream.
     Returns:
      val - value of variable in file.
     Arguments:
      FID - stream to read from.
      name - name of string.
        

    See also
    writeStringToFID, readIntFromFID, readBoolFromFID, readDoubleFromFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    lineStr = getline(FID);
    tokens = tokenise(lineStr, '=');
    if(length(tokens)~=2 | ~strcmp(tokens{1}, string))
      error('Incorrect file format.')
    end
    val = tokens{2};

def val = readVersionFromFID(FID):

    """Read version number from an FID.
    
    Description:
    
    val = readVersionFromFID(FID) reads version number from a stream.
     Returns:
      val - value of variable in file.
     Arguments:
      FID - stream to read from.
        

    See also
    writeVersionToFID, readDoubleFromFID, readStringFromFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    val = str2num(readStringFromFID(FID, 'version'));

def [area, rocPointX, rocPointY] = rocCurve(outputs, labels):

    """Draw ROC curve and return labels.
    
    Description:
    
    [area, rocPointX, rocPointY] = rocCurve(outputs, labels) draws an
     ROC curve and returns the area under the ROC curve as well as the
     points plotted.
     Returns:
      area - teh area under the ROC curve.
      rocPointX - the x points of the ROC curve.
      rocPointY - the y points of the ROC curve.
     Arguments:
      outputs - the outputs from the model (e.g. probabilities of
       labels).
      labels - the true labels associated with the outputs.
        

    See also
    trapz


    Copyright (c) 2004 Neil D. Lawrence
    
    """
        
    
    [sOutputs, index] = sort(outputs);
    sLabels = labels(index);
    
    for i = length(sOutputs):-1:1;
      % False positives
      rocPointX(length(sOutputs)-i+1) = sum(sLabels(i:end)==-1)/sum(labels==-1);
      rocPointY(length(sOutputs)-i+1) = sum(sLabels(i:end)==1)/sum(labels==1);
    end
    
    if nargin < 3
      plot(rocPointX, rocPointY);
    end
    area = trapz(rocPointX, rocPointY);
def [ax, data] = scatterPlot(X, YLbls, markerSize):

    """2-D scatter plot of labelled points.
    
    Description:
    
    [ax, data] = scatterPlot(X, lbls, markerSize) plots a 2-D scatter
     plot of labelled points.
     Returns:
      ax - the axis handle of the plot.
      data - the handles of the data points.
     Arguments:
      X - x points to plot.
      lbls - the labels to plot.
      markerSize - the size of the markers to use.
        

    See also
    plot, getSymbols


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    if isempty(YLbls)
      symbol = [];
    else
      symbol = getSymbols(size(YLbls,2));
    end
    
    % Create the plot for the data
    clf
    ax = axes('position', [0.05 0.05 0.9 0.9]);
    hold on
    
    data = twoDPlot(X, YLbls, symbol);
    for i = 1:2
      minX = min(X(:, i));
      maxX = max(X(:, i));
      xSpan = maxX - minX;
      xLim{i}(1) = minX - 0.1*xSpan;
      xLim{i}(2) = maxX + 0.1*xSpan;
    end
    if nargin>2
      for i = 1:length(data)
        set(data(i), 'markersize', markerSize);
      end
    end
    set(ax, 'xLim', xLim{1});
    set(ax, 'yLim', xLim{2});
    
    set(ax, 'fontname', 'arial');
    set(ax, 'fontsize', 20);
    
def returnVal = twoDPlot(X, label, symbol):
    
    % GPLVMTWODPLOT Helper function for plotting the labels in 2-D.
    
    returnVal = [];
    
    if ~isempty(label)
      for i = 1:size(X, 1)
        labelNo = find(label(i, :));
        try 
          returnVal = [returnVal; plot(X(i, 1), X(i, 2), symbol{labelNo})];
        catch
          if strcmp(lasterr, 'Index exceeds matrix dimensions.')
    	error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
          end
        end
      end
    else
      returnVal = plot(X(:, 1), X(:, 2), 'rx');
    end
def y = sigmoid(x):

    """The sigmoid function
    
    Description:
    def y = sigmoid(x):
%
    """
        
    y = ones(size(x))./(1+exp(-x));
def D = sparseDiag(d):

    """Create a diagonal matrix that is sparse from a vector.
    
    Description:
    
    D = sparseDiag(d) creates a diagonal matrix that is sparse from a
     vector.
     Returns:
      D - the sparse diagonal matrix containing the vector as its
       diagonal.
     Arguments:
      d - the diagonal vector from which the sparse diagonal matrix is
       formed.
        

    See also
    diag, spdiags


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    if length(size(d)) ~=2
      error('Input must be a vector.');
    end
    if size(d, 1) ~= 1 & size(d, 2) ~=1
      error('Input must be a vector.');
    end
    
    D = spdiags(d, 0, length(d), length(d));
    % % Can be made more efficient.
    % n = length(d);
    % D = spalloc(n, n, n);
    % for i = 1:n
    %   D(i, i) = d(i);
    % end

def x = stack(X):

    """Return column stacked vector of given matrix.
    
    Description:
    
    x = stack(X) returns a column stacked vector of a given matrix.
     This is useful if you wish to stack a column vector from a matrix
     returned by another function (i.e. you can't apply the colon
     operator directly).
     Returns:
      x - stacked column vector from the matrix.
     Arguments:
      X - the matrix to be stacked.
        

    See also
    colon


    Copyright (c) 2004 Neil D. Lawrence
    
    """
        
    
    x = X(:);
def str = stringSigFigs(num, sf):

    """Convert number to a string with a number of significant digits.
    
    Description:
    
    str = stringSigFigs(number, sf) converts a given number to a
     string, but provides a given number of significant figures.
     Returns:
      str - the string with the number to the given number of
       significant figures.
     Arguments:
      number - the number that requires conversion.
      sf - the number of significant figures required in the conversion.
        

    See also
    num2str


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    val = chop(num, sf);
    str = num2str(val);
    ind = 1;
    while str(ind) == '0' | str(ind) == '.'
      ind = ind+1;
    end
    count = 0;
    while(ind <= length(str))
      if str(ind) ~= '.'
        count = count + 1;
      end
      ind = ind +1;
    end
    while count < 3
      str = [str '0'];
      count = count + 1;
    end

def parts = stringSplit(string, separator):

    """Return separate parts of a string.
    
    Description:
    
    parts = stringSplit(string, separator) return separate parts of a
     string split by a given separator.
     Returns:
      parts - cell array of the string split into parts.
     Arguments:
      string - the string to be split into parts.
      separator - the character used to split the string (default, ',').
        

    See also
    tokenise


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    if nargin < 2
      separator = ',';
    end
    
    parts = tokenise(string, separator);
def [columnNames, data] = tableRead(fileName, separator):

    """Read in data which has column titles in the first line and separated values in each other line.
    
    Description:
    
    [columnNames, data] = tableRead(fileName, separator) reads in data
     from a file that has column titles in the first line and separated
     values in every other line.
     Returns:
      columnNames - the names of the columns taken from the first line.
      data - the data, taken from the remaining lines.
     Arguments:
      fileName - file name in which the data is stored.
      separator - separator between the columns (default ',').
        

    See also
    fopen, fgetl


    Copyright (c) 2004 Neil D. Lawrence
    
    """
        
    if nargin < 2
      separator = ',';
    end
    
    fid = fopen(fileName);
    i = 0;
    while 1
      i = i + 1;
      lin=fgetl(fid);
      if ~ischar(lin), break, end
    end
    numLines = i;
    
    fid = fopen(fileName);
    lin = fgetl(fid);
    columnNames = stringSplit(lin, separator);
    numCol = length(columnNames);
    i = 0;
    data = zeros(numLines - 1, numCol);
    while 1
      i = i+1;
      lin=fgetl(fid);
      if ~ischar(lin), break, end
      split = stringSplit(lin, separator);
      if length(split) ~= numCol
        error(['Error at line ' num2str(i) ' of file ' fileName ': wrong ' ...
                            'number of columns'])
      end
      for j = 1:length(split);
        data(i, j) = num2str(split{j});
      end
    end
    fclose(fid);
def tokens = tokenise(string, delim):

    """Split a string into separate tokens.
    
    Description:
    
    tokens = tokenise(string, delim) takes a string and splits it into
     separate tokens.
     Returns:
      tokens - a cell array of tokens split from the string.
     Arguments:
      string - the string to be split into parts.
      delim - the delimiter to use in splitting the string (default is '
       ').
        

    See also
    % SEEALSO stringSplit


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    if nargin < 2
      delim = ' ';
    end
    tokpos = find(string==delim);
    len=length(tokpos);
    
    if len==0
      tokens{1} = string;
    elseif len==1
      tokens{1} = string(1:tokpos(1)-1);
      tokens{2} = string(tokpos(1)+1:end);
    elseif len>1
      tokens{1} = string(1:tokpos(1)-1);
      for(i=2:len)
        tokens{i} = string(tokpos(i-1)+1:tokpos(i)-1);
      end
      tokens{len+1} = string(tokpos(len)+1:end);
    end
def t = traceProduct(A, B):

    """Returns the trace of the product of two matrices.
    
    Description:
    
    t = traceProduct(A, B) returns the trace of the product of two
     matrices, tr(A*B).
     Returns:
      t - the trace of the product.
     Arguments:
      A - the first matrix in the product.
      B - the second matrix in the product.
        

    See also
    trace


    Copyright (c) 2004 Neil D. Lawrence
    
    """
        
    t = sum(sum(A.*B'));
def tree = treeFindChildren(tree):

    """Given a tree that lists only parents, add children.
    
    Description:
    
    tree = treeFindChildren(tree) takes a tree structure which lists
     the children of each node and computes the parents for each node
     and places them in.
     Returns:
      tree - a tree that lists children and parents.
     Arguments:
      tree - the tree that lists only children.
        

    See also
    treeFindParents


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    for i = 1:length(tree)
      for j = 1:length(tree(i).parent)
        if tree(i).parent(j)
          tree(tree(i).parent(j)).children ...
              = [tree(tree(i).parent(j)).children i];
        end
      end
    end
    

def ind = treeFindLeaves(tree):

    """Return indices of all leaf nodes in a tree structure.
    
    Description:
    
    ind = treeFindLeaves(tree) returns indices of all leaf nodes in an
     tree array.
     Returns:
      ind - indices of leaf nodes.
     Arguments:
      tree - tree for which leaf nodes are being sought.
        

    See also
    treeFindParents, treeFindChildren


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    ind = [];
    for i = 1:length(tree)
      if isempty(tree(i).children)
        ind = [ind i];
      end
    end
def tree = treeFindParents(tree):

    """Given a tree that lists only children, add parents.
    
    Description:
    
    tree = treeFindParents(tree) takes a tree structure which lists
     the parents of each node and computes the children for each node
     and places them in.
     Returns:
      tree - a tree that lists parents and children.
     Arguments:
      tree - the tree that lists only parents.
        

    See also
    treeFindChildren


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    for i = 1:length(tree)
      for j = 1:length(tree(i).children)
        if tree(i).children(j)
          tree(tree(i).children(j)).parent ...
              = [tree(tree(i).children(j)).parent i];
        end
      end
    end
    

def ind = treeFindRoots(tree):

    """Return indices of all root nodes in a tree structure.
    
    Description:
    
    ind = treeFindRoots(tree) returns indices of all root nodes in an
     tree array.
     Returns:
      ind - indices of root nodes.
     Arguments:
      tree - tree for which root nodes are being sought.
        

    See also
    treeFindParents, treeFindChildren, treeFindLeaves


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    ind = [];
    for i = 1:length(tree)
      if isempty(tree(i).parent)
        ind = [ind i];
      end
    end
def [widths, maxDepth, nodePositions] = treeGetWidths(tree):

    """give width of each level of tree.
    
    Description:
    
    [widths, maxDepth, nodePositions] = treeGetWidths(tree) gives the
     width of a tree at each level of the hierarchy and the maximum
     depth of a tree as well as the node stored at a give depth and
     breadth.
     Returns:
      widths - stores the width at each depth level.
      maxDepth - the maximum depth of the tree.
      nodePositions - stores the nodeIndex present at depth i and
       breadth j.
     Arguments:
      tree - the tree for which the dimensions are required.
        

    See also
    treeFindParents, treeFindChildren


    Copyright (c) 2006 Andrew J. Moore
    
    """
        
    maxDepth = 0;
    %widths stores the width at each depth level. Since the max depth isn't yet
    %known, allocate enough memory to widths for the worst case scenario where
    %each node has only 1 child and the depth is equal to the number of nodes.
    widths = zeros(length(tree), 1);
    %nodePositions 
    nodePositions = zeros(length(tree), length(tree));
    rootInd = treeFindRoots(tree);
    indAlreadySeen = [];
    for i = 1:length(rootInd)
      traverseTree(rootInd(i), 1);
    end
    widths = widths(1:maxDepth, :); %trim off excess rows
    
      function traverseTree(nodeIndex, depthLevel)
        if ~any(indAlreadySeen == nodeIndex)
          widths(depthLevel) = widths(depthLevel) + 1;
        end
        indAlreadySeen = [indAlreadySeen nodeIndex];
    
        nodePositions(depthLevel, widths(depthLevel)) = nodeIndex;
        if length(tree(nodeIndex).children) > 0
          for i=1:length(tree(nodeIndex).children)
            traverseTree(tree(nodeIndex).children(i), depthLevel + 1);
          end
        else
          if depthLevel > maxDepth
            maxDepth = depthLevel;
          end
        end
      end
    
    end
def tree = treeSwapNode(tree, i, j):

    """Swap two nodes in the tree structure array.
    
    Description:
    
    tree = treeSwapNode(tree, i, j) swaps the location of two nodes in
     a tree structure array.
     Returns:
      tree - the tree structure with the two node locations swapped.
     Arguments:
      tree - the tree for which two nodes are to be swapped.
      i - the index of the first node to be swapped.
      j - the index of the second node to be swapped.
        

    See also
    treeFindParents, treeFindChildren


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    storeNodeI = tree(i);
    storeNodeJ = tree(j);
    tree(j) = storeNodeI;
    tree(i) = storeNodeJ;
    for k = 1:length(tree)
      tree(k).children(find(tree(k).children==i)) = -1;
      tree(k).children(find(tree(k).children==j)) = i;
      tree(k).children(find(tree(k).children==-1)) = j;
      tree(k).parent(find(tree(k).parent==i)) = -1;
      tree(k).parent(find(tree(k).parent==j)) = i;
      tree(k).parent(find(tree(k).parent==-1)) = j;
    end

def writeBoolToFID(FID, name, val):

    """Writes a boolean to an FID.
    
    Description:
    
    writeBoolToFID(FID, name, val) writes a boolean to a stream.
     Arguments:
      FID - stream to write to.
      name - name of boolean.
      val - value of variable to place in file.
        

    See also
    readBoolFromFID, writeStringToFID, writeBoolToFID, writeIntToFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    if val
      writeStringToFID(FID, name, '1');
    else
      writeStringToFID(FID, name, '0');
    end
      

def writeDoubleToFID(FID, name, val):

    """Writes a double to an FID.
    
    Description:
    
    writeDoubleToFID(FID, name, val) writes a double to a stream.
     Arguments:
      FID - stream to write to.
      name - name of double.
      val - value of variable to place in file.
        

    See also
    readDoubleFromFID, writeStringToFID, writeBoolToFID, writeIntToFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    writeStringToFID(FID, name, num2str(val));

def writeIntToFID(FID, name, val):

    """Writes an integer to an FID.
    
    Description:
    
    writeIntToFID(FID, name, val) writes an integer to a stream.
     Arguments:
      FID - stream to write to.
      name - name of int.
      val - value of variable to place in file.
        

    See also
    readIntFromFID, writeStringToFID, writeBoolToFID, writeDoubleToFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    writeStringToFID(FID, name, num2str(val));

def writeStringToFID(FID, name, val):

    """Writes a string to an FID.
    
    Description:
    
    writeStringToFID(FID, name, val) writes an string from a stream.
     Arguments:
      FID - stream to write from.
      name - name of string.
      val - value of variable to place in file.
        

    See also
    readStringFromFID, writeIntToFID, writeBoolToFID, writeDoubleToFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    fprintf(FID, [name '=' val '\n']);

def writeVersionToFID(FID, val):

    """Writes a version to an FID.
    
    Description:
    
    writeVersionToFID(FID, val) writes a version from a stream.
     Arguments:
      FID - stream to write to.
      val - value of version to place in file.
        

    See also
    readVersionFromFID, writeStringToFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    writeStringToFID(FID, 'version', num2str(val));
      

def z = xlogy(x, y):

    """z = x*log(y) returns zero if x=y=0
    
    Description:
    
    z = xlogy(x, y) computes the function x*log(y) but taking account
     for the fact that the answer is zero if x=y=0.
     Returns:
      z - returns z where z = x*log(y).
     Arguments:
      x - the x argument in x*log(y).
      y - the y argument in x*log(y).
    
    y = xlogy(x) computes the function x*log(x), taking account for
     the fact that the answer is zero if x=0.
     Returns:
      y - returns y where y = x*log(x).
     Arguments:
      x - the argument in x*log(x).


    Copyright (c) 2001, 2006 Neil D. Lawrence
    
    """
        
    % If there is only one input argument x is assumed to equal y.
    if nargin == 1
      y=x;
    end
    if any(any(x==0))
      if ~issparse(x)
        z = zeros(size(x));
      else
        z = spalloc(size(x, 1), size(x, 2), sum(sum(x~=0)));
      end
      indx = find(x);
      z(indx)= x(indx).*log(y(indx));
    else
      z= x.*log(y);
    end
def zeroAxes(axesHandle, tickRatio, fontSize, fontName):

    """A function to move the axes crossing point to the origin.
    
    Description:
    
    zeroAxes(axesHandle, tickRatio, fontSize, fontName) moves the
     crossing point of the axes to the origin.
     Arguments:
      axesHandle - the handle of the axes to zero.
      tickRatio - the ratio of the axis length to the tick length
       (default 0.025).
      fontSize - the font size for the axes (default 14).
      fontName - the name of the font for the axes (default 'times').
        

    See also
    plot


    Copyright (c) 2005 Neil D. Lawrence
    
    """
        
    if nargin < 4
      fontName = [];
      if nargin < 3
        fontSize = [];
        if nargin < 2
          tickRatio = [];
          if nargin < 1
            axesHandle = [];
          end
        end
      end
    end
    if isempty(fontName)
      fontName = 'times';
    end
    if isempty(fontSize)
      fontSize = 14;
    end
    if isempty(tickRatio)
      tickRatio = 0.025;
    end
    if isempty(axesHandle)
      axesHandle = [];
    end
    axis off
    xlim = get(axesHandle, 'xlim');
    ylim = get(axesHandle, 'ylim');
    xTickLength = (xlim(2) - xlim(1))*tickRatio;
    yTickLength = (ylim(2) - ylim(1))*tickRatio;
    
    xtick = get(axesHandle, 'xtick');
    ytick = get(axesHandle, 'ytick');
    xaxYpos = 0;
    if ylim(1) > 0
      xaxYpos = ylim(1);
    elseif ylim(2) < 0
      xaxYpos = ylim(2);
    end
    line([xlim(1) xlim(2)], [xaxYpos xaxYpos])
    for i = 1:length(xtick)
      if xtick(i) == 0
        if any(ytick == 0)
          continue
        end
      end
      line([xtick(i) xtick(i)], xaxYpos+[0 -yTickLength]);
      textHandle = text(xtick(i), xaxYpos-4*yTickLength, num2str(xtick(i)));
      set(textHandle, 'fontsize', fontSize, 'fontname', fontName, ...
    		  'horizontalalignment', 'center');
    end
    
    yaxXpos = 0;
    if xlim(1) > 0
      yaxXpos = xlim(1);
      set(axesHandle, 'xlim', [xlim(1) - xTickLength xlim(2)]);
    elseif xlim(2) < 0
      yaxXpos = xlim(2);
      set(axesHandle, 'xlim', [xlim(1) + xTickLength xlim(2)]);
    end
    line([yaxXpos  yaxXpos], [ylim(1) ylim(2)])
    for i = 1:length(ytick)
      if ytick(i) == 0
        if any(xtick ==  0)
          t = 0:pi/24:2*pi;
          xCirc = xTickLength*sin(t);
          yCirc = yTickLength*cos(t);
          line(xCirc', yCirc');
          continue
        end
      end
      line(yaxXpos-[0 xTickLength], [ytick(i) ytick(i)]);
      textHandle = text(yaxXpos - 1.1*xTickLength, ytick(i),  num2str(ytick(i)));
      set(textHandle, 'fontsize', fontSize, 'fontname', fontName, ...
    		  'horizontalalignment', 'right');
    
    end
    
    
    
    

