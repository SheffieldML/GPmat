% GPLVMOIL Model the oil data with a 2-D GPLVM.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

load 3class

Y = DataTrn;
lbls = DataTrnLbls;

% Centre the data
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
X = Y*u(:, 1:2)*diag(1./sqrt(v(1:2)));


returnVal = [];

symbol{1} = 'r+';
symbol{2} = 'bo';
symbol{3} = 'mx';
figure, hold on
for i = 1:size(X, 1)
  labelNo = find(lbls(i, :));
  plot(X(i, 1), X(i, 2), symbol{labelNo})
end


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
gplvmvisualise(X, Y, invK, theta, lbls, meanData, activeSet, 'vectorVisualise', ...
	       'vectorModify');
save gplvmOil.mat X theta activeSet