function [X, y, XTest, yTest] = ivmLoadData(dataset)

% IVMLOADDATA Load a dataset.

% IVM

randn('seed', 1e5)
rand('seed', 1e5)
switch dataset
 case 'usps'
  load ../data/usps_train
  X = ALL_DATA;
  range =  min(ALL_T):max(ALL_T);
  for i = 1:length(range)
    y(:, i) = (ALL_T == range(i))*2 - 1;
  end
  if nargout > 2
    load ../data/usps_test
    XTest = ALL_DATA;
    range =  min(ALL_T):max(ALL_T);
    for i = 1:length(range)
      yTest(:, i) = (ALL_T == range(i))*2 - 1;
    end
  end

 case 'regressionOne'
  try
    load regressionOneData.mat 
  catch
    
    [void, errid] = lasterr;
    if strcmp(errid, '')
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
    
    save('regressionOneData.mat', 'numIn', 'N', 'X', 'y')
  end
  
 case 'regressionTwo'
  try
    load regressionTwoData.mat 
  catch
    
    [void, errid] = lasterr;
    if strcmp(errid, '')
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
      kern.inputScales = [0.999 .2];
      
      K = kernCompute(kern, X);
      y = real(gaussSamp(K, 1)') + randn(N, 1)*0.01;
    end
    
    save('regressionTwoData.mat', 'numIn', 'N', 'X', 'y')
  end

 case 'classificationTwo'
   
  try
    load classificationTwoData.mat 
  catch
    
    noFile = 0;
    verString = version;
    
    if str2double(verString(1:3)) > 6.1
      [void, errid] = lasterr;
      if strcmp(errid, '')
        noFile = 1;
      end
    else
      errMsg = lasterr;
      if findstr(errMsg, 'file does not exist')
        noFile = 1;
      end
    end
    if noFile
      numIn = 2;
      N = 500;
      
      X = zeros(N, numIn);
      X = rand(N, numIn);
      
      kern = kernCreate(X, 'rbf');
      kern.variance = 100;
      kern.inverseWidth = 10;
      
      K = kernCompute(kern, X);
      u = real(gaussSamp(K, 1)');
      
      p = cumGaussian(u);
      y = 2*(rand(size(u))>p)-1;
    end
    
    save('classificationTwoData.mat', 'numIn', 'N', 'X', 'u', 'y')
  end
 
 case 'orderedOne'
  dataPerCat = 30;
  spacing = 3;
  
  % Generate a toy data-set of (linear) ordered categories.
  X = [randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]; ...
       randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
       randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
       randn(dataPerCat, 2); ...
       randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
       randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
       randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]];
  y = [zeros(dataPerCat, 1); ...
       repmat(1, dataPerCat, 1); repmat(2, dataPerCat, 1); ...
       repmat(3, dataPerCat, 1); repmat(4, dataPerCat, 1); ...
       repmat(5, dataPerCat, 1); repmat(6, dataPerCat, 1)];

 case 'orderedTwo'
  dataPerCat = 30;
  spacing = 3;
  
  % Generate a toy data-set of (linear) ordered categories.
  thetaR = [randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
            randn(dataPerCat, 2); ...
            randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]];
  thetaR(:, 1) = thetaR(:, 1);
  thetaR(:, 2) = thetaR(:, 2) + 15;
  X = [sin(thetaR(:, 1)).*thetaR(:, 2) cos(thetaR(:, 1)).*thetaR(:, 2)];
  y = [zeros(dataPerCat, 1); ...
       repmat(1, dataPerCat, 1); repmat(2, dataPerCat, 1); ...
       repmat(3, dataPerCat, 1); repmat(4, dataPerCat, 1); ...
       repmat(5, dataPerCat, 1); repmat(6, dataPerCat, 1)];

end
