function [X, y, XTest, yTest] = mappingLoadData(dataset, seedVal)

% MAPPINGLOADDATA Load a regression or classification dataset.


if nargin < 2
  seedVal = 1e5;
end
randn('seed', seedVal)
rand('seed', seedVal)
XTest = [];
yTest = [];
switch dataset
 %/~
 case 'pumadynSeeger'
  data = load('Dataset.data');
  ind = randperm(size(data, 1));
  indTr = ind(1:7168);
  indTe = ind(7169:end);
  X = data(indTr, 1:end-1);
  y = data(indTr, end);
  XTest = data(indTe, 1:end-1);
  yTest = data(indTe, end);
  Xscale = sqrt(var(X));
  for j = 1:length(Xscale);
    X(:, j) = X(:, j)/Xscale(j);
    XTest(:, j) = XTest(:, j)/Xscale(j);
  end
  augX = [X ones(size(X, 1), 1)];
  w = inv(augX'*augX)*augX'*y;
  y = y -augX*w;
  augTestX = [XTest ones(size(XTest, 1), 1)];
  yTest = yTest - augTestX*w;
  yscale = sqrt(var(y));
  y = y/yscale;
  yTest = yTest/yscale;
 case 'pumadyn'

  % Data is variance 1, no need to normalise.
  data = load('Dataset.data');
  ind = randperm(size(data, 1));
  indTr = ind(1:7168);
  indTe = ind(7169:end);
  X = data(indTr, 1:end-1);
  y = data(indTr, end);
  XTest = data(indTe, 1:end-1);
  yTest = data(indTe, end);
  %~/
 case 'usps'
  load usps_train
  X = ALL_DATA;
  range =  min(ALL_T):max(ALL_T);
  for i = 1:length(range)
    y(:, i) = (ALL_T == range(i))*2 - 1;
  end
  if nargout > 2
    load usps_test
    XTest = ALL_DATA;
    range =  min(ALL_T):max(ALL_T);
    for i = 1:length(range)
      yTest(:, i) = (ALL_T == range(i))*2 - 1;
    end
  end
  
 case {'usps0', 'usps1', 'usps2', 'usps3', 'usps4', 'usps5', 'usps6', 'usps7', 'usps8', 'usps9'}
  digitNo = str2num(dataset(end));
  load usps_train
  X = ALL_DATA;
  range =  min(ALL_T):max(ALL_T);
  for i = 1:length(range)
    y(:, i) = (ALL_T == range(i))*2 - 1;
  end
  if nargout > 2
    load usps_test
    XTest = ALL_DATA;
    range =  min(ALL_T):max(ALL_T);
    for i = 1:length(range)
      yTest(:, i) = (ALL_T == range(i))*2 - 1;
    end
  end
  y = y(:, digitNo+1);
  yTest = yTest(:, digitNo+1);
  
 case 'regressionOne'
  try
    load regressionOneData.mat 
  catch
    
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
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

      save('regressionOneData.mat', 'numIn', 'N', 'X', 'y')
    else
      error(lasterr);
    end
    
  end
  
 case 'regressionTwo'
  try
    load regressionTwoData.mat 
  catch
    
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
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
    
      save('regressionTwoData.mat', 'numIn', 'N', 'X', 'y')
    else
      error(lasterr)
    end
        
  end

 case 'regressionThree'
  try
    load regressionThreeData.mat 
  catch
    
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
      numIn = 2;
      N = 500;
      X = zeros(N, numIn);
      X(1:floor(N/2), :) = ...
          randn(floor(N/2), numIn)*.5+1;
      X(floor(N/2)+1:end, :) = ...
          randn(ceil(N/2), numIn)*.5-1;
      kern = kernCreate(X, 'lin');
      kern.variance = 1;
      
      K = kernCompute(kern, X);
      y = real(gaussSamp(K, 1)') + randn(N, 1)*0.01;
      save('regressionThreeData.mat', 'numIn', 'N', 'X', 'y')
    else 
      error(lasterr);
    end
  end
  
 case 'regressionFour'
  try
    load regressionFourData.mat 
  catch
    
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
      numIn = 2;
      N = 500;
      X = zeros(N, numIn);
      X(1:floor(N/2), :) = ...
          randn(floor(N/2), numIn)*.5+1;
      X(floor(N/2)+1:end, :) = ...
          randn(ceil(N/2), numIn)*.5-1;
      kern = kernCreate(X, 'mlp');
      kern.variance = 1;
      kern.weightVariance = 1;
      kern.biasVariance = 1;
      K = kernCompute(kern, X);
      y = real(gaussSamp(K, 1)') + randn(N, 1)*0.01;
      save('regressionFourData.mat', 'numIn', 'N', 'X', 'y')
    else 
      error(lasterr);
    end
  end
 case 'classificationOne'
  try
    load classificationOneData.mat
  catch
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
      X = [randn(100,2)-[zeros(100, 1) 6*ones(100, 1)]; randn(100,2)+[zeros(100, 1) 6*ones(100, 1)]; randn(100, 2)];
      y = [ones(200, 1); -ones(100, 1)];
      save('classificationOneData.mat', 'X', 'y')
    else
      error(lasterr);  
    end
  end
 case 'classificationTwo'
   
  try
    load classificationTwoData.mat 
  catch
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
      numIn = 2;
      N = 500;
      
      X = zeros(N, numIn);
      X = rand(N, numIn);
      
      kern = kernCreate(X, 'rbf');
      kern.variance = 10;
      kern.inverseWidth = 10;
      
      K = kernCompute(kern, X);
      u = real(gaussSamp(K, 1)');
      
      p = cumGaussian(u);
      y = 2*(rand(size(u))>p)-1;
      save('classificationTwoData.mat', 'numIn', 'N', 'X', 'u', 'y')
    else
      error(lasterr);  
    end
    
  end

 case 'classificationThree'
   
  try
    load classificationThreeData.mat 
  catch
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
      numIn = 2;
      N = 500;
      
      X = zeros(N, numIn);
      X = rand(N, numIn);
      
      kern = kernCreate(X, 'rbf');
      kern.variance = 10;
      kern.inverseWidth = 10;
      
      K = kernCompute(kern, X);
      u = real(gaussSamp(K, 1)');
      a = 3;
      pMinus = cumGaussian(u-a/2);
      pPlus = cumGaussian(-u-a/2);
      p =rand(size(u));
      indMinus = find(p<pMinus);
      indPlus = find(p>pMinus & p<pMinus+pPlus);
      indNone = find(p>pMinus+pPlus);
      y = zeros(N, 1);
      y(indPlus) = 1;
      y(indMinus) = -1;
      y(indNone, :) = [];
      X(indNone, :) = [];
      save('classificationThreeData.mat', 'numIn', 'N', 'X', 'u', 'y')
    else
      error(lasterr);  
    end
    
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
