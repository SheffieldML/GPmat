% GENERATEREGRESSIONDATA Tries to load a sampled data set otherwise generates it.

% IVM

try
  load regressionData.mat 
catch
  
  [void, errid] = lasterr;
  if strcmp(errid, '')
    randn('seed', 1e5)
    rand('seed', 1e5)
    numIn = 4;
    N = 2000;
    trueTheta = [10 1 0.01 0 1 .5 1 .5 1];
    trueTheta = thetaConstrain(trueTheta);
    trueLntheta = log(trueTheta);
    
    X = zeros(N, numIn);
    X(1:floor(N/2), :) = ...
      randn(floor(N/2), numIn)*.25+1;
    X(floor(N/2)+1:end, :) = ...
	randn(ceil(N/2), numIn)*.25-1;
  
    K = kernel(X, trueLntheta, 'ARD');
    y = gaussSamp(K, 1)';
  end

  save('regressionData.mat', 'numIn', 'N', 'X', 'y')
end
