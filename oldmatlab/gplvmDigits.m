% GPLVMDIGITS Model the digits data with a 2-D GPLVM.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% load data
load c:\datasets\usps\matlab\16x16\usps_train.mat

% Extract 600 of digits 0 to 4
[ALL_T, sortIndices] = sort(ALL_T);
ALL_DATA = ALL_DATA(sortIndices(:), :);
Y = [];
lbls = [];
numEachDigit = 600;
for digit = 0:4;
  firstDigit = min(find(ALL_T==digit));
  Y = [Y; ALL_DATA(firstDigit:firstDigit+numEachDigit-1, :)];
  lbl = zeros(1, 5);
  lbl(digit+1) = 1;
  lbls = [lbls; repmat(lbl, numEachDigit, 1)];
end

% Don't Centre the data
meanData = zeros(1, size(Y, 2)); %mean(Y);
%Y = Y  - repmat(meanData, size(Y, 1), 1);

% Extract data dimensions and set IVM active set size
numData = size(Y, 1);
dataDim = size(Y, 2);
numActive = 100;
extIters = 15;

% Initialise X with PCA
[v, u] = pca(Y);
v(find(v<0))=0;
X = Y*u(:, 1:2)*diag(1./sqrt(v(1:2)));

% Initialise theta
theta(1) = 1;
theta(2) = 1;
theta(3) = 1;

% Options for optimisation in latent space
options = foptions;
options(1) = 0;
options(9) = 0;
options(14) = 100;

% options for kernel optimisation
optionsKernel = foptions;
optionsKernel(1) = 0;
optionsKernel(9) = 0;
optionsKernel(14) = 100;

% Fit the GP latent variable model
[X, theta, activeSet] = gplvmfit(X, Y, theta, numActive, optionsKernel, ...
		      options, extIters)
  
% compute the kernel from results
[K, invK] = computeKernel(X(activeSet, :), theta);

% Visualise the results
gplvmvisualise(X, Y, invK, theta, lbls, meanData, activeSet, 'imageVisualise', ...
	       'imageModify', [16 16]);
%save gplvmDigits.mat X theta activeSet