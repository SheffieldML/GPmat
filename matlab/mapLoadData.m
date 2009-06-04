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
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

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

 case 'silhouette'
  % Ankur Agarwal and Bill Trigg's silhoutte data.
  load([baseDir dirSep 'mocap' dirSep 'ankur' dirSep 'ankurDataPoseSilhouette']);
  inMean = mean(Y);
  inScales = sqrt(var(Y));
  X = Y - repmat(inMean, size(Y, 1), 1);
  X = X./repmat(inScales, size(Y, 1), 1);
  
  XTest = Y_test - repmat(inMean, size(Y_test, 1), 1);
  XTest = XTest./repmat(inScales, size(Y_test, 1), 1);
  y = Z;
  yTest = Z_test;
  
  
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
  %/~
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
  %~/
 case 'ggToy'
  try
    load([baseDir 'toyMultigp1D.mat']);
  catch
    randn('seed', 1e5); % This was taken to make the outputs equal to the ones of the NIPS paper
    rand('seed', 1e5);  % This was taken to make the outputs equal to the ones of the NIPS paper
    nout =4;
    prec = [50 50 300 200];
    sigma2_y = [1 1 5 5];
    media = zeros(1,4);
    prec_u = 100;
    sigma2_u = 1;
    N = 500;
    N2 = 100;
    x = linspace(-1,1,N)';
    x2 = linspace(-1,1,N2)';
    Kyy = cell(nout);
    ggKern1.inputDimension = size(x,2);
    ggKern1.precision_u = prec_u;
    ggKern1.sigma2_u = sigma2_u;
    ggKern2.inputDimension = size(x,2);
    for i = 1:nout,
      ggKern1.precision_y = prec(i);
      ggKern1.sigma2_y = sigma2_y(i);
      ggKern1.translation = media(i) ;
      for j = 1:nout,
        ggKern2.precision_y = prec(j);
        ggKern2.sigma2_y = sigma2_y(j);
        ggKern2.translation = media(j);
        Kyy{i,j} = ggXggKernCompute(ggKern1, ggKern2, x);
      end
    end
    Kyu = cell(nout, 1);
    ggKern.inputDimension = size(x,2);
    gaussianKern.precision_u = prec_u;
    gaussianKern.sigma2_u = sigma2_u;
    for i = 1:nout,
      ggKern.precision_y = prec(i);
      ggKern.translation = media(i) ;
      ggKern.sigma2_y = sigma2_y(i);
      Kyu{i,j} = ggXgaussianKernCompute(ggKern, gaussianKern, x, x2);
    end
    Kuu = cell(1);
    gaussianKern.precision_u = prec_u;
    gaussianKern.sigma2_u = sigma2_u;
    Kuu{1} = gaussianKernCompute(gaussianKern, x2);
    K = [cell2mat(Kuu) cell2mat(Kyu)';cell2mat(Kyu) cell2mat(Kyy)];
    yu = gsamp(zeros(size(K,1),1), K, 1);
    y = yu(size(x2,1)+1:end) ;
    Y = reshape(y,size(x,1),nout);
    for k=1:nout,
      Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
    end
    missingData = cell(nout,1);
    missingData{1} = 101:160;
    missingData{4} = 21:90;
    ntrainx =200;
    randn('seed', 1e4); % This was taken to make the outputs equal to the ones of the NIPS paper
    rand('seed', 1e4);  % This was taken to make the outputs equal to the ones of the NIPS paper
    maxl = length(x);
    X = cell(1,nout);
    y = cell(1,nout);
    yTest = cell(1,nout);
    indx = randperm(maxl);
    pindx = sort(indx(1:ntrainx));
    for k =1:nout,
      X{k} = x(pindx,:);
      X{k}(missingData{k},:)= [];
      y{k} = Y(pindx,k);
      y{k}(missingData{k}) = [];
      yTest{k} = Y(indx(ntrainx+1:end),k);
    end
    XTest = x(indx(ntrainx+1:end),:)';
    save([baseDir 'toyMultigp1D.mat'], 'X', 'y', 'XTest', 'yTest')
  end
    case 'ggToyCombined'
        try
            load([baseDir 'ggToyCombined.mat']);
        catch
            noise = 1;
            nout =5;
            nin = 2;
            media = zeros(1,5);
            prec = [50 10 300 100 30];
            sigma2_y = [1 1 5 5 2];
            sens_y_noise = [1 2 2.5 10 1];
            prec_u = 100;
            sigma2_u = 1;
            N = 500;
            N2 = 100;
            x = linspace(-1,1,N)';
            x2 = linspace(-1,1,N2)';
            Kyy = cell(nout);
            ggKern1.inputDimension = size(x,2);
            ggKern1.precision_u = prec_u;
            ggKern1.sigma2_u = sigma2_u;
            ggKern1Noise.inputDimension =size(x,2);
            ggKern1Noise.sigma2Noise = noise;
            ggKern2.inputDimension = size(x,2);
            ggKern2Noise.inputDimension = size(x,2);
            for i = 1:nout,
                ggKern1.precision_y = prec(i);
                ggKern1.sigma2_y = sigma2_y(i);
                ggKern1.translation = media(i) ;
                ggKern1Noise.precisionG = prec(i);
                ggKern1Noise.variance = sens_y_noise(i);
                for j = 1:nout,
                    ggKern2.precision_y = prec(j);
                    ggKern2.sigma2_y = sigma2_y(j);
                    ggKern2.translation = media(j);
                    ggKern2Noise.precisionG = prec(j);
                    ggKern2Noise.variance = sens_y_noise(j);
                    Kyy{i,j} = ggXggKernCompute(ggKern1, ggKern2, x) + ggwhiteXggwhiteKernCompute(ggKern1Noise, ...
                        ggKern2Noise, x);
                end
            end
            Kyu = cell(nout, 2);
            ggKern.inputDimension = size(x,2);
            gaussianKern.precision_u = prec_u;
            gaussianKern.sigma2_u = sigma2_u;
            ggKernNoise.inputDimension = size(x,2);
            for i = 1:nout,
                ggKern.precision_y = prec(i);
                ggKern.translation = media(i) ;
                ggKern.sigma2_y = sigma2_y(i);
                Kyu{i,1} = ggXgaussianKernCompute(ggKern, gaussianKern, x, x2);
            end
            for i = 1:nout,
                ggKernNoise.precisionT = 2*prec(i);
                ggKernNoise.sigma2Noise = noise;
                Kyu{i,2} = sens_y_noise(i)*gaussianwhiteKernCompute(ggKernNoise, x, x2);
            end
            Kuu = cell(2,1);
            gaussianKern.precision_u = prec_u;
            gaussianKern.sigma2_u = sigma2_u;
            Kuu{1} = gaussianKernCompute(gaussianKern, x2);
            Kuu{2} = noise*eye(length(x2));
            startVal = 1;
            endVal = 0;
            KuuMat = zeros(2*length(x2));
            for k = 1:nin,
                endVal = endVal + length(x2);
                KuuMat(startVal:endVal,startVal:endVal) = Kuu{k};
                startVal = endVal + 1;
            end
            K = [KuuMat cell2mat(Kyu)';cell2mat(Kyu) cell2mat(Kyy)];
            yu = gsamp(zeros(size(K,1),1), K, 1);
            u = yu(1:2*size(x2,1));
            y = yu(2*size(x2,1)+1:end);
            U = reshape(u,size(x2,1),nin);
            Y = reshape(y,size(x,1),nout);
            for k=1:nout,
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            ntrainx =200;
            maxl = length(x);
            X = cell(1,nout+nin);
            y = cell(1,nout+nin);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = x(pindx,:);
                y{k} = Y(pindx,k);
                yTest{k} = Y(indx(ntrainx+1:end),k);
            end
            XTest = x(indx(ntrainx+1:end),:)';
            % Append the latent functions
            y{nout+1} = U(:,1);
            y{nout+2} = U(:,2);
            X{nout+1} = x2;
            X{nout+2} = x2;
            save([baseDir 'ggToyCombined.mat'], 'X', 'y', 'XTest', 'yTest')
        end
    case 'ggToyCombinedMissing'
        try
            load([baseDir 'ggToyCombinedMissing.mat']);
        catch
            [X, y, XTest, yTest] = mapLoadData('ggToyCombined');
            nout = 5;         
            missingData = cell(nout,1);
            missingData{1} = 37:90;
            %missingData{3} = 108:148;
            missingData{5} = 147:184;
            for k =1:nout,
                X{k}(missingData{k},:)= [];
                y{k}(missingData{k}) = [];               
            end
            save([baseDir 'ggToyCombinedMissing.mat'], 'X', 'y', 'XTest', 'yTest')
        end
    case 'ggwhiteToy'
        try
            load([baseDir 'ggwhiteToy.mat']);
        catch
            noise = 1;
            nout =5;
            nin = 1;           
            prec = [50 10 300 100 30];
            sens_y_noise = [1 2 2.2 1.5 1];
            N = 500;
            N2 = 100;
            x = linspace(-1,1,N)';
            x2 = linspace(-1,1,N2)';
            Kyy = cell(nout);
            ggKern1Noise.inputDimension =size(x,2);
            ggKern1Noise.sigma2Noise = noise;
            ggKern2Noise.inputDimension = size(x,2);
            for i = 1:nout,
                ggKern1Noise.precisionG = prec(i);
                ggKern1Noise.variance = sens_y_noise(i);
                for j = 1:nout,
                    ggKern2Noise.precisionG = prec(j);
                    ggKern2Noise.variance = sens_y_noise(j);
                    Kyy{i,j} = ggwhiteXggwhiteKernCompute(ggKern1Noise, ggKern2Noise, x);
                end
            end
            Kyu = cell(nout, nin);
            ggKernNoise.inputDimension = size(x,2);
            ggKernNoise.variance = noise;
            for i = 1:nout,                
                ggKern1Noise.variance = sens_y_noise(i);
                ggKern1Noise.precisionG = prec(i);                
                Kyu{i,1} = ggwhiteXwhiteKernCompute(ggKern1Noise, ggKernNoise, x, x2);
            end
            Kuu = cell(nin,1);            
            Kuu{1} = noise*eye(length(x2));
            K = [cell2mat(Kuu) cell2mat(Kyu)';cell2mat(Kyu) cell2mat(Kyy)];
            yu = gsamp(zeros(size(K,1),1), K, 1);
            u = yu(1:size(x2,1));
            y = yu(size(x2,1)+1:end);
            U = reshape(u,size(x2,1),nin);
            Y = reshape(y,size(x,1),nout);
            for k=1:nout,
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            ntrainx =200;
            maxl = length(x);
            X = cell(1,nout+nin);
            y = cell(1,nout+nin);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = x(pindx,:);
                y{k} = Y(pindx,k);
                yTest{k} = Y(indx(ntrainx+1:end),k);
            end
            XTest = x(indx(ntrainx+1:end),:)';
            % Append the latent functions
            y{nout+1} = U(:,1);
            X{nout+1} = x2;            
            save([baseDir 'ggwhiteToy.mat'], 'X', 'y', 'XTest', 'yTest')
        end
        
    case 'ggwhiteToyMissing'
        try
            load([baseDir 'ggwhiteToyMissing.mat']);
        catch
            [X, y, XTest, yTest] = mapLoadData('ggwhiteToy');
            nout = 5;         
            missingData = cell(nout,1);
            missingData{1} = 46:108;
            missingData{4} = 147:184;            
            for k =1:nout,
                X{k}(missingData{k},:)= [];
                y{k}(missingData{k}) = [];               
            end
            save([baseDir 'ggwhiteToyMissing.mat'], 'X', 'y', 'XTest', 'yTest')
        end
    case 'ggwhiteToyHighDim'
        try
            load([baseDir 'ggwhiteToyHighDim.mat']);
        catch
            noise = 1;
            nout =5;
            nin = 1;           
            inputDim = 2;
            prec = [50 10 5 100 30];
            sens_y_noise = [1 2 2.2 1.5 1];
            N = 20;
            N2 = 20;
            x = linspace(-1,1,N)';
            x2 = linspace(-1,1,N2)';
            [Xout1, Xout2 ] = ndgrid(x);
            [Xlf1, Xlf2 ] = ndgrid(x2);
            Xout = [Xout1(:) Xout2(:)]; 
            Xlf = [Xlf1(:)  Xlf2(:)];            
            Kyy = cell(nout);
            ggKern1Noise.inputDimension = inputDim;
            ggKern1Noise.sigma2Noise = noise;
            ggKern2Noise.inputDimension = inputDim;
            for i = 1:nout,
                ggKern1Noise.precisionG = prec(i)*ones(inputDim,1);
                ggKern1Noise.variance = sens_y_noise(i);
                for j = 1:nout,
                    ggKern2Noise.precisionG = prec(j)*ones(inputDim,1);
                    ggKern2Noise.variance = sens_y_noise(j);
                    Kyy{i,j} = ggwhiteXggwhiteKernCompute(ggKern1Noise, ggKern2Noise, Xout);
                end
            end
            Kyu = cell(nout, nin);
            ggKernNoise.inputDimension = inputDim;
            ggKernNoise.variance = noise;
            for i = 1:nout,                
                ggKern1Noise.variance = sens_y_noise(i);
                ggKern1Noise.precisionG = prec(i)*ones(inputDim,1);                
                Kyu{i,1} = ggwhiteXwhiteKernCompute(ggKern1Noise, ggKernNoise, Xout, Xlf);
            end
            Kuu = cell(nin,1);            
            Kuu{1} = noise*eye(size(Xlf,1));
            K = [cell2mat(Kuu) cell2mat(Kyu)';cell2mat(Kyu) cell2mat(Kyy)];
            yu = gsamp(zeros(size(K,1),1), K, 1);
            u = yu(1:size(Xlf,1));
            y = yu(size(Xlf,1)+1:end);
            U = reshape(u,size(Xlf,1),nin);
            Y = reshape(y,size(Xout,1),nout);
            for k=1:nout,
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            ntrainx =200;
            maxl = size(Xout,1);
            X = cell(1,nout+nin);
            y = cell(1,nout+nin);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = Xout(pindx,:);
                y{k} = Y(pindx,k);
                yTest{k} = Y(indx(ntrainx+1:end),k);
            end
            XTest = Xout(indx(ntrainx+1:end),:)';
            % Append the latent functions
            y{nout+1} = U(:,1);
            X{nout+1} = Xlf;            
            save([baseDir 'ggwhiteToyHighDim.mat'], 'X', 'y', 'XTest', 'yTest')
        end
        
    case 'ggwhiteToyMissingHighDim'
        try
            load([baseDir 'ggwhiteToyMissingHighDim.mat']);
        catch
            [X, y, XTest, yTest] = mapLoadData('ggwhiteToyHighDim');
            nout = 5;         
            missingData = cell(nout,1);
            missingData{1} = 46:108;
            missingData{4} = 147:184;            
            for k =1:nout,
                X{k}(missingData{k},:)= [];
                y{k}(missingData{k}) = [];               
            end
            save([baseDir 'ggwhiteToyMissingHighDim.mat'], 'X', 'y', 'XTest', 'yTest')
        end        
    case 'simToyCombined'
        try
            load([baseDir 'simToyCombined.mat']);
        catch
            noise = 1;
            nout =5;
            nin = 2;
            decay = [2 10 5 0.1 7.5];
            variance = [1.5 1.5 4 10 2];
            sensitivity = [1 0.3 1.2 0.2 0.5];
            inverseWidth = 100;
            N = 200;
            N2 = 50;
            x = linspace(0,1,N)';
            x2 = linspace(0,1,N2)';
            Kyy = cell(nout);
            simKern1.inputDimension = size(x,2);
            simKern1 = simKernParamInit(simKern1);
            simKern1.inverseWidth = inverseWidth;
            simwhiteKern1.inputDimension =size(x,2);
            simwhiteKern1 = simwhiteKernParamInit(simwhiteKern1);
            simwhiteKern1.variance = noise;
            simKern2.inputDimension = size(x,2);
            simKern2 = simKernParamInit(simKern2);
            simKern2.inverseWidth = inverseWidth;
            simwhiteKern2.inputDimension = size(x,2);
            simwhiteKern2 = simwhiteKernParamInit(simwhiteKern2);
            simwhiteKern2.variance = noise;
            for i = 1:nout,
                simKern1.variance = variance(i);
                simKern1.decay = decay(i);
                simwhiteKern1.decay = decay(i);
                simwhiteKern1.sensitivity = sensitivity(i);
                for j = 1:nout,
                    simKern2.variance = variance(j);
                    simKern2.decay = decay(j);                    
                    simwhiteKern2.decay = decay(j);
                    simwhiteKern2.sensitivity = sensitivity(j);
                    Kyy{i,j} = simXsimKernCompute(simKern1, simKern2, x) + ...
                        simwhiteXsimwhiteKernCompute(simwhiteKern1, simwhiteKern2, x);
                end
            end
            Kyu = cell(nout, 2);
            rbfKern.inputDimension = size(x,2);
            rbfKern = rbfKernParamInit(rbfKern);
            rbfKern.inverseWidth = inverseWidth;
            for i = 1:nout,
                simKern1.variance = variance(i);
                simKern1.decay = decay(i);                
                Kyu{i,1} = simXrbfKernCompute(simKern1, rbfKern, x, x2);
            end
            whiteKern.inputDimension = size(x,2);            
            whiteKern = whiteKernParamInit(whiteKern);
            whiteKern.variance = noise;
            for i = 1:nout,
                simwhiteKern1.decay = decay(i);
                simwhiteKern1.sensitivity = sensitivity(i);
                Kyu{i,2} = simwhiteXwhiteKernCompute(simwhiteKern1, whiteKern, x, x2);
            end
            Kuu = cell(2,1);
            Kuu{1} = rbfKernCompute(rbfKern, x2);
            Kuu{2} = noise*eye(length(x2));
            startVal = 1;
            endVal = 0;
            KuuMat = zeros(2*length(x2));
            for k = 1:nin,
                endVal = endVal + length(x2);
                KuuMat(startVal:endVal,startVal:endVal) = Kuu{k};
                startVal = endVal + 1;
            end
            K = [KuuMat cell2mat(Kyu)';cell2mat(Kyu) cell2mat(Kyy)];
            yu = real(gsamp(zeros(size(K,1),1), K, 1));
            u = yu(1:2*size(x2,1));
            y = yu(2*size(x2,1)+1:end);
            U = reshape(u,size(x2,1),nin);
            Y = reshape(y,size(x,1),nout);
            for k=1:nout,
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            ntrainx =100;
            maxl = length(x);
            X = cell(1,nout+nin);
            y = cell(1,nout+nin);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = x(pindx,:);
                y{k} = Y(pindx,k);
                yTest{k} = Y(indx(ntrainx+1:end),k);
            end
            XTest = x(indx(ntrainx+1:end),:)';
            % Append the latent functions
            y{nout+1} = U(:,1);
            y{nout+2} = U(:,2);
            X{nout+1} = x2;
            X{nout+2} = x2;
            save([baseDir 'simToyCombined.mat'], 'X', 'y', 'XTest', 'yTest')
        end
    case 'ggwhiteToyHighDimBatch'
        try
            load([baseDir 'ggwhiteToyHighDimBatch.mat']);
        catch
            kernName = 'ggwhite';            
            inverseWidth = [50 10 5   100 30 200 400 20 1 40];
            sensitivity =  [1  2  2.2 1.5 1  3   4   2  1    10];            
            noisePerOutput = 1e-2;
            nlf = 1;
            nout = 1;
            d = nlf + nout;
            inputDim = 5;
            % Sample the inputs              
            X1 = gsamp(zeros(inputDim,1), eye(inputDim), 200);
%            X1 = linspace(-1,1, 200)';
            kernType{1} = multigpKernComposer(kernName, d, nlf, 'ftc', 1);
            kernType{2} = multigpKernComposer('white',  d, nlf, 'ftc', 1);
            kern = kernCreate(X1,  {'cmpnd', kernType{:}});
            kern.comp{1}.comp{1}.variance = 1;
            for k = 1:nout,
                kern.comp{1}.comp{1+k}.precisionG = inverseWidth(k)*ones(inputDim,1);
                kern.comp{1}.comp{1+k}.variance = sensitivity(k);
                kern.comp{2}.comp{1+k}.variance = noisePerOutput;
            end
            K = kernCompute(kern, X1);
            yu = gsamp(zeros(size(K,1),1), K, 1);
            u = yu(1:size(X1,1));
            y = yu(size(X1,1)+1:end);            
            Kout = K(1+size(X1,1):end,1+size(X1,1):end);            
            [invKout, U, jitter] = pdinv(Kout);
            if any(jitter>1e-4)
                fprintf('Warning: Added jitter of %2.4f\n', jitter)
            end
            logDetKout = logdet(Kout, U);
            dim = size(y, 2);
            ll = -dim*log(2*pi) -logDetKout - y*invKout*y';
            ll = ll*0.5;
            U = reshape(u,size(X1,1),nlf);
            Y = reshape(y,size(X1,1),nout);           
            X = cell(1, nout);
            y = cell(1, nout);            
            for k =1:nout,
                X{k} = X1;
                y{k} = Y(:,k);
            end
            % Append the latent functions
            y{nout+1} = U(:,1);
            X{nout+1} = X1;
            XTest = [];
            yTest = [];
            save([baseDir 'ggwhiteToyHighDimBatch.mat'], 'X', 'y', 'XTest', 'yTest')
        end
                            
    case 'compilerData'
         data = load([baseDir 'data_compiler_org.mat']);
         X = cell(1,length(data.X));
         y = cell(1,length(data.X));
%          for k = 1:length(data.X),
%              X{k} = zscore(data.X{k});
%              y{k} = zscore(data.Y{k});
%          end   
         X = data.X';
         y= data.Y';
         XTest = [];
         yTest = [];
         
    case 'schoolData'
        data = load([baseDir 'schoolData.mat']);
        X = cell(1,length(data.X));
        y = cell(1,length(data.X));
%         for k = 1:length(data.X),
%             X{k} = zscore(data.X{k});
%             y{k} = zscore(data.y{k});
%         end
        X = data.X';
        y= data.y';
        XTest = [];
        yTest = [];

    case 'juraDataCd'
  try
    load([baseDir 'juraDataCd.mat']);
  catch           
    fidP = fopen([baseDir 'prediction.dat'],'r');
    if fidP ==-1
      error('The file prediction.dat does not exist in this directory');
    end
    fidV = fopen([baseDir 'validation.dat'],'r');
    if fidV ==-1
      error('The file validation.dat does not exist in this directory');
    end
    xLoc = [1 2];
    primaryLoc = 5; % Column location of Cd in file prediction.dat and file validation.dat
    secondaryLoc = [9 11]; % Column location of Ni and Zn 
    fseek(fidP, 61,-1);
    A = textscan(fidP, '%f%f%f%f%f%f%f%f%f%f%f');
    predValues = cell2mat(A);
    fseek(fidV, 61,-1);
    A = textscan(fidV, '%f%f%f%f%f%f%f%f%f%f%f');
    valValues = cell2mat(A);
    X = cell(1, 1+length(secondaryLoc));
    y = cell(1, 1+length(secondaryLoc));
    XTest =  cell(1, 1+length(secondaryLoc));
    yTest =  cell(1, 1+length(secondaryLoc));
    X{1} = predValues(:, xLoc);
    y{1} = predValues(:,primaryLoc);
    for i=1:length(secondaryLoc);
      X{1+i} = predValues(:, xLoc);
      y{1+i} = predValues(:, secondaryLoc(i));
    end
    % Append validation values to training data
    for i=1:length(secondaryLoc);
      X{1+i} = [X{1+i}; valValues(:, xLoc)];
      y{1+i} = [y{1+i}; valValues(:, secondaryLoc(i))];
    end
    % Form the testting sets
    XTest{1} = valValues(:, xLoc);
    yTest{1} = valValues(:, primaryLoc);
    for i=1:length(secondaryLoc);
      XTest{1+i} = valValues(:, xLoc);
      yTest{1+i} = valValues(:, secondaryLoc(i));
    end
    fclose(fidP);
    fclose(fidV);
    save([baseDir 'juraDataCd.mat'], 'X', 'y', 'XTest', 'yTest');
    
  end
        %/~
 case 'juraDataCo'
  data = load([baseDir 'data_jura_Co']);
  X = data.Xtrain;
  y = data.Ytrain;
  XTest = data.Xtest;
  yTest = data.Ytest;
  %~/
 case 'juraDataCu'
  try
    load([baseDir 'juraDataCu.mat']);
  catch
    fidP = fopen([baseDir 'prediction.dat'],'r');
    if fidP ==-1
      error('The file prediction.dat does not exist in this directory');
    end
    fidV = fopen([baseDir 'validation.dat'],'r');
    if fidV ==-1
      error('The file validation.dat does not exist in this directory');
    end
    xLoc = [1 2];
    primaryLoc = 8; % Column location of Cu in file prediction.dat and file validation.dat
    secondaryLoc = [10 9 11]; % Column location of Ni and Zn
    fseek(fidP, 61,-1);
    A = textscan(fidP, '%f%f%f%f%f%f%f%f%f%f%f');
    predValues = cell2mat(A);
    fseek(fidV, 61,-1);
    A = textscan(fidV, '%f%f%f%f%f%f%f%f%f%f%f');
    valValues = cell2mat(A);
    X = cell(1, 1+length(secondaryLoc));
    y = cell(1, 1+length(secondaryLoc));
    XTest =  cell(1, 1+length(secondaryLoc));
    yTest =  cell(1, 1+length(secondaryLoc));
    X{1} = predValues(:, xLoc);
    y{1} = predValues(:,primaryLoc);
    for i=1:length(secondaryLoc);
      X{1+i} = predValues(:, xLoc);
      y{1+i} = predValues(:, secondaryLoc(i));
    end
    % Append validation values to training data
    for i=1:length(secondaryLoc);
      X{1+i} = [X{1+i}; valValues(:, xLoc)];
      y{1+i} = [y{1+i}; valValues(:, secondaryLoc(i))];
    end
    % Form the testting sets
    XTest{1} = valValues(:, xLoc);
    yTest{1} = valValues(:, primaryLoc);
    for i=1:length(secondaryLoc);
      XTest{1+i} = valValues(:, xLoc);
      yTest{1+i} = valValues(:, secondaryLoc(i));
    end
    fclose(fidP);
    fclose(fidV);
    save([baseDir 'juraDataCu.mat'], 'X', 'y', 'XTest', 'yTest');
  end
  %/~
 case 'juraDataPb'
  data = load([baseDir 'data_jura_Pb']);
  X = data.Xtrain;
  y = data.Ytrain;
  XTest = data.Xtest;
  yTest = data.Ytest;
  
 case 'demp53_5genes'
  data = load([baseDir 'dataBarencoOption_0_Genes_5']);
  X = data.data{1}.Xtrain;
  y = [cell2mat(data.data{1}.Ytrain); cell2mat(data.data{2}.Ytrain); cell2mat(data.data{3}.Ytrain)];
  XTest = X;
  yTest = y;
  
 case 'demp53_50genes'
  data = load([baseDir 'dataBarencoOption_0_Genes_50']);
  X = data.data{1}.Xtrain;
  y = [cell2mat(data.data{1}.Ytrain); cell2mat(data.data{2}.Ytrain); cell2mat(data.data{3}.Ytrain)];
  XTest = X;
  yTest = y;
  
 case 'demp53_20genes'
  data = load([baseDir 'dataBarencoOption_0_Genes_20']);
  X = data.data{1}.Xtrain;
  y = [cell2mat(data.data{1}.Ytrain); cell2mat(data.data{2}.Ytrain); cell2mat(data.data{3}.Ytrain)];
  XTest = X;
  yTest = y;
  %~/
 otherwise
  error('Unknown data set requested.')
  
end

if semiSup % Test if data is for semi-supervised learning.
  indUnlabelled = find(rand(size(y, 1), 1)>labProb);
  y(indUnlabelled, :) = NaN;
end