% GENERATEREGRESSIONDATA Tries to load a sampled data set otherwise generates it.


try
  load regressionData.mat 
catch
  
  [void, errid] = lasterr;
  if strcmp(errid, '')
    randn('seed', 1e5)
    rand('seed', 1e5)
    numIn = 2;
    N = 500;
    X = zeros(N, numIn);
    X(1:floor(N/2), :) = ...
      randn(floor(N/2), numIn)*.5+1;
    X(floor(N/2)+1:end, :) = ...
	randn(ceil(N/2), numIn)*.5-1;
    kern = kernCreate(X, 'rbfard');
    kern.variance = 1;
    kern.inverseWidth = 20;
    kern.inputScales = [0 0.999];
    
    K = kernCompute(kern, X);
    y = real(gaussSamp(K, 1)') + randn(N, 1)*0.01;
  end

  save('regressionData.mat', 'numIn', 'N', 'X', 'y')
end
