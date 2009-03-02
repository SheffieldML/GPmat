function [X, y, XTest, yTest] = mapLoadData(dataset, seedVal)

% MAPLOADDATA Load a mapping model dataset (e.g. classification, regression).
% FORMAT
% DESC loads a data set for a mapping modelling problem
% (e.g. classification or regression).
% ARG dataset : the name of the data set to be loaded. Currently
% the possible names are 'usps', 'usps.' where . is a digit between
% 0 and 9, 'spgp1d', 'twoclusters', 'regressionOne',
% 'regressionTwo', 'regressionThree', 'regressionFour',
% 'classificationTwo', 'classificationThree', 'orderedOne',
% 'orderedTwo', 'ionosphere'.
% ARG seedVal : a seed value for generating the data set (default
% is 1e5). Note that in many cases a generated data set will be
% saved and loaded from disc the next time the function is
% called. In these cases this value, if changed, won't have an effect.
% RETURN X : the training input data loaded in.
% RETURN y : the training target data loaded in.
% RETURN XTest : the test input data loaded in. If no test set is
% available it is empty.
% RETURN yTest : a test target data.
%
% SEEALSO : lvmLoadData, datasetsDirectory
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006


% DATASETS

if nargin < 2
  seedVal = 1e5;
end
randn('seed', seedVal)
rand('seed', seedVal)
XTest = [];
yTest = [];

baseDir = datasetsDirectory;
dirSep = filesep;
semiSup = false;

if length(dataset)>5 && strcmp(dataset(1:6), 'gunnar')
  % Data set is one of Gunnar Raetsch's
  ind = find(dataset==':');
  dataSetName = dataset(ind(1)+1:ind(2)-1);
  dataSetNum = dataset(ind(2)+1:end);
  filebase = [baseDir dirSep 'gunnar' dirSep dataSetName dirSep dataSetName];
  X=load([filebase '_train_data_' num2str(dataSetNum) '.asc']);
  y=load([filebase '_train_labels_' num2str(dataSetNum) '.asc']);
  XTest=load([filebase '_test_data_' num2str(dataSetNum) '.asc']);
  yTest=load([filebase '_test_labels_' num2str(dataSetNum) '.asc']);
  return
elseif length(dataset)>3 && strcmp(dataset(1:4), 'semi')
  % Data set is semi-supervised learning.
  ind = find(dataset==':');
  labProb = str2num(dataset(ind(2)+1:end));
  % Extract dataset part
  dataset = dataset(ind(1)+1:ind(2)-1);
  semiSup = true;
end
  
switch dataset
  
 case 'cedar69'
  % Data for ICML 2001 paper on noisy KFD.
  digOne = 6;
  digTwo = 9;

  directory = [baseDir '\cedar.cd\matlab\16x16\'];
  load([directory 'digit' num2str(digOne)]);
  X = digitData;
  numDataOne = size(digitData, 1);
  y = zeros(size(digitData, 1), 1);
  mydigitData = [];
  digitData = [];
  load([directory 'digit' num2str(digTwo)]);
  X = [X; digitData];
  X = double(X);
  X = -(X-128)/128; 
  numDataTwo = size(digitData, 1);
  y = [y; ones(size(digitData, 1), 1)];
  
  if nargin>2
    % Load test data
    load([directory 'bsdigit' num2str(digOne)]);
    XTest = digitData;
    numDataOne = size(digitData, 1);
    yTest = zeros(size(digitData, 1), 1);
    digitData = [];
    load([directory 'bsdigit' num2str(digTwo)]);
    XTest = [XTest; digitData];
    XTest = double(XTest);
    XTest = -(XTest-128)/128; 
    yTest = [yTest; ones(size(digitData, 1), 1)];
  end
 
 case 'sky'
  load([baseDir 'skydatatrain.mat']);
  X = double(x)/255;
  y = double(t)*2 -1;
  load([baseDir 'skydatavalid.mat']);
  XTest = double(x)/255;
  yTest = double(t)*2 -1;
  
 case 'unlabelledOne'
  numDataPart = 100;
  [X, y] = generateCrescentData(numDataPart);
  
 case 'classificationOne'
  try
    load([baseDir 'classificationOneData'])
  catch
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
      X = [randn(100,2)-[zeros(100, 1) 6*ones(100, 1)]; randn(100,2)+[zeros(100, 1) 6*ones(100, 1)]; randn(100, 2)];
      y = [ones(200, 1); -ones(100, 1)];
      save([baseDir 'classificationOneData.mat'], 'X', 'y')
    else
      error(lasterr);  
    end
  end
 case 'usps'
  load([baseDir 'usps_train']);
  X = ALL_DATA;
  range =  min(ALL_T):max(ALL_T);
  for i = 1:length(range)
    y(:, i) = (ALL_T == range(i))*2 - 1;
  end
  if nargout > 2
    load([baseDir 'usps_test']);
    XTest = ALL_DATA;
    range =  min(ALL_T):max(ALL_T);
    for i = 1:length(range)
      yTest(:, i) = (ALL_T == range(i))*2 - 1;
    end
  end
  
 case {'usps0', 'usps1', 'usps2', 'usps3', 'usps4', 'usps5', 'usps6', 'usps7', 'usps8', 'usps9'}
  digitNo = str2num(dataset(end));
  load([baseDir 'usps_train'])
  X = ALL_DATA;
  range =  min(ALL_T):max(ALL_T);
  for i = 1:length(range)
    y(:, i) = (ALL_T == range(i))*2 - 1;
  end
  if nargout > 2
    load([baseDir 'usps_test']);
    XTest = ALL_DATA;
    range =  min(ALL_T):max(ALL_T);
    for i = 1:length(range)
      yTest(:, i) = (ALL_T == range(i))*2 - 1;
    end
  end
  y = y(:, digitNo+1);
  yTest = yTest(:, digitNo+1);
 case 'threeFive'
   load([baseDir 'usps_train'])
   X = ALL_DATA;
   y = ALL_T;
   load([baseDir 'usps_test'])
   XTest = ALL_DATA;
   yTest = ALL_T;
   classTrue = 3;
   for i = [0 1 2 4 6 7 8 9];
     index = find(y == i);
     X(index, :) = [];
     y(index, :) = [];
     index = find(yTest == i);
     XTest(index, :) = [];
     yTest(index, :) = [];
   end
   y = (y == classTrue)*2 - 1;
   yTest = (yTest == classTrue)*2 - 1;
 case 'fourNine'
  load([baseDir 'usps_train'])
  X = ALL_DATA;
  y = ALL_T;
  load([baseDir 'usps_test'])
  XTest = ALL_DATA;
  yTest = ALL_T;
  classTrue = 4;
  for i = [0 1 2 3 5 6 7 8];
    index = find(y == i);
    X(index, :) = [];
    y(index, :) = [];
    index = find(yTest == i);
    XTest(index, :) = [];
    yTest(index, :) = [];
  end
  y = (y == classTrue)*2 - 1;
  yTest = (yTest == classTrue)*2 - 1;
  
 case 'thorsten'
  [y, X] = svmlread([baseDir 'example2/train_transduction.dat']);
  [yTest, XTest] = svmlread([baseDir 'example2/test.dat']);

 case 'spgp1d'
  try 
    load([baseDir 'spgp1DData.mat'])
  catch
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
      numIn = 1;
      N = 500;
      X = 2*rand(N, numIn)-1;
      kern = kernCreate(X, 'rbf');
      kern.variance = 1;
      kern.inverseWidth = 20;      
      K = kernCompute(kern, X);
      y = real(gaussSamp(K, 1)') + randn(N, 1)*0.1;

      save([baseDir 'spgp1DData.mat'], 'numIn', 'N', 'X', 'y')
    else
      error(lasterr);
    end
  end
 case 'twoclusters'
  try 
    load([baseDir 'twoclusters.mat'])
  catch
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
      numIn = 1;
      N = 200;
      X1 = rand(N/2, numIn)-1.5;
      X2 = rand(N/2, numIn)+.5;
      X = [X1; X2];
      kern = kernCreate(X, 'rbf');
      kern.variance = 1;
      kern.inverseWidth = 100;      
      K = kernCompute(kern, X);
      y = real(gaussSamp(K, 1)') + randn(N, 1)*0.1;

      save([baseDir 'twoclusters.mat'], 'numIn', 'N', 'X', 'y')
    else
      error(lasterr);
    end
  end
 case 'regressionOne'
  try
    load([baseDir 'regressionOneData.mat'])
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

      save([baseDir 'regressionOneData.mat'], 'numIn', 'N', 'X', 'y')
    else
      error(lasterr);
    end
    
  end
  
 case 'regressionTwo'
  try
    load([baseDir 'regressionTwoData.mat'])
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
    
      save([baseDir 'regressionTwoData.mat'], 'numIn', 'N', 'X', 'y')
    else
      error(lasterr)
    end
        
  end

 case 'regressionThree'
  try
    load([baseDir 'regressionThreeData.mat'])
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
      save([baseDir 'regressionThreeData.mat'], 'numIn', 'N', 'X', 'y')
    else 
      error(lasterr);
    end
  end

 case 'regressionFour'
  try
    load([baseDir 'regressionFourData.mat'])
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
      save([baseDir 'regressionFourData.mat'], 'numIn', 'N', 'X', 'y')
    else 
      error(lasterr);
    end
  end

 case 'classificationTwo'
   
  try
    load([baseDir 'classificationTwoData.mat'])
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
      save([baseDir 'classificationTwoData.mat'], 'numIn', 'N', 'X', 'u', 'y')
    else
      error(lasterr);  
    end
    
  end

 case 'classificationThree'
   
  try
    load([baseDir 'classificationThreeData.mat'])
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
      indPlus = find(p>pMinus && p<pMinus+pPlus);
      indNone = find(p>pMinus+pPlus);
      y = zeros(N, 1);
      y(indPlus) = 1;
      y(indMinus) = -1;
      y(indNone, :) = [];
      X(indNone, :) = [];
      save([baseDir 'classificationThreeData.mat'], 'numIn', 'N', 'X', 'u', 'y')
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


 case 'ionosphere'
  X = zeros(351, 34);
  fid = fopen([baseDir 'ionosphere.data'], 'r');
  lin = getline(fid);
  i = 0;
  while(lin ~= -1)
    i = i+1;
    elements = tokenise(lin, ',');
    for j = 1:length(elements)-1
      X(i, j) = str2num(elements{j});
    end
    switch(elements{end})
     case 'g'
      y(i, 1) = 1;
     case 'b'
      y(i, 1) = -1;
    end
    lin = getline(fid);
  end
  ind = randperm(351);
  trainInd = ind(1:200);
  testInd = ind(201:end);
  XTest = X(testInd, :);
  yTest = y(testInd, :);
  X = X(trainInd, :);
  y = y(trainInd, :);

  %/~
 case 'pumadynSeeger'
  data = load([baseDir 'Dataset.data']);
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
  data = load([baseDir 'Dataset.data']);
  ind = randperm(size(data, 1));
  indTr = ind(1:7168);
  indTe = ind(7169:end);
  X = data(indTr, 1:end-1);
  y = data(indTr, end);
  XTest = data(indTe, 1:end-1);
  yTest = data(indTe, end);
 %~/
 
    case 'lfmOde'      
      data = load([baseDir 'datasetODE31_150_5Ind0']);
      X = data.Xtrain{1};
      y = cell2mat(data.Ytrain');
      XTest = data.Xtest{1};
      yTest = cell2mat(data.Ytest'); 
      
    case 'lfmOdeInd'      
      data = load([baseDir 'datasetODE31_100_5Ind1']);
      X = data.Xtrain{1};
      y = cell2mat(data.Ytrain');
      XTest = data.Xtest{1};
      yTest = cell2mat(data.Ytest'); 
        
    case 'ggToy'
        data = load([baseDir 'sample141']);
        X = data.X';
        y = data.Y';
        XTest = data.Xtest';
        yTest = data.Ytest';
    
    case 'juraDataCd'    
        data = load([baseDir 'data_jura_Cd']);
        X = data.Xtrain;
        y = data.Ytrain;
        XTest = data.Xtest;
        yTest = data.Ytest;
    
    case 'juraDataCo'    
        data = load([baseDir 'data_jura_Co']);
        X = data.Xtrain;
        y = data.Ytrain;
        XTest = data.Xtest;
        yTest = data.Ytest;
        
    case 'juraDataCu'    
        data = load([baseDir 'data_jura_Cu']);
        X = data.Xtrain;
        y = data.Ytrain;
        XTest = data.Xtest;
        yTest = data.Ytest;
        
    case 'juraDataPb'    
        data = load([baseDir 'data_jura_Pb']);
        X = data.Xtrain;
        y = data.Ytrain;
        XTest = data.Xtest;
        yTest = data.Ytest;    
              
 otherwise
  error('Unknown data set requested.')

end

if semiSup % Test if data is for semi-supervised learning.
  indUnlabelled = find(rand(size(y, 1), 1)>labProb);
  y(indUnlabelled, :) = NaN;
end