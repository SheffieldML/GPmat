% GPLVMBRENDAN1D Model the face data with a 1-D GPLVM.
% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% load data
load c:\datasets\faces\brendan\frey_rawface.mat

Y = double(ff)';

% Find mean brendan
meanData = mean(Y);
Y = Y  - repmat(meanData, size(Y, 1), 1);

% Extract data dimensions and set IVM active set size
numData = size(Y, 1);
dataDim = size(Y, 2);
numActive = 100;
extIters = 15;

% Initialise X with PCA

[v, u] = pca(Y);
v(find(v<0))=0;
X = Y*u(:, 1)*diag(1./sqrt(v(1)));

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
gplvmvisualise1D(X, Y, invK, theta, [], meanData, activeSet, 'imageVisualise', ...
		 'imageModify', [20 28]);

% save gplvmBrendan1D.mat X theta activeSet