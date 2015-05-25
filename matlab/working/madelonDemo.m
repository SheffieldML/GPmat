% MADELONDEMO Try the IVM with ARD kernel on the MADELON data-set.

randn('seed', 1e5)
rand('seed', 1e5)

model.kernelType = 'ARD';
prior = 1;
display = 1;
rawData = dlmread('c:\datasets\challenge\MADELON\madelon_train.data', ' ');
rawData = rawData(:, 1:end-1);
meanData = mean(rawData);
rawData = rawData - repmat(meanData, size(rawData, 1), 1);
sdData = sqrt(var(rawData));
X = rawData./repmat(sdData, size(rawData, 1), 1);

kurt =kurtosis(X);
[void, order] = sort(abs(kurt));
model.X = X(:, order(end-50:end));
model.y = dlmread('c:\datasets\challenge\MADELON\madelon_train.labels', ' ');

numData = size(model.X, 1);
switch model.kernelType
 case 'regular'
  model.lntheta = [zeros(1, 5)];
 case 'ARD'
  model.lntheta = [zeros(1, 5) zeros(1, size(model.X, 2))];
end

%model.lntheta = factor*[ones(1, 5) ...
%		    log([0 0 0 1 0 0 0 0 0 ...
%		    0 0 1 0 0 0 0 0 1 ...
%		    0 0 0 0 0 0 0 0 0 ...
%		    0 0 0 0 0 0 0 0 0 ...
%		    0 0 0 0 1 1 1 1 1 ...
%		    1 1 1 0 0 1]+1e-16)];


% Number of inclusions
model.d = 600;

nClass1 = sum(model.y==1);
nClass2 = sum(model.y== -1);

% Set the bias parameter
model.bias = invCummGaussian(nClass1/(nClass2+nClass1));

% selection criteria
selectionCriteria = 'entropy';
for j = 1:4
  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseKernel(model, prior, display);
end
disp(model.lntheta)
model.d = 600
model = ivmOptimiseIVM(model, display);

xTest = dlmread('c:\datasets\challenge\MADELON\madelon_valid.data', ' ');
yTest = dlmread('c:\datasets\challenge\MADELON\madelon_valid.labels', ' ');
yPred = ivmfwd(xTest, model);

testError = 1-sum(yPred==yTest)/size(yTest, 1);