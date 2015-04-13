% GPLVRESULTSMOIL Load and visualise the results for the oil algorithm.

load 3class

Y = DataTrn;
lbls = DataTrnLbls;

% Centre the data
meanData = mean(Y);
Y = Y  - repmat(meanData, size(Y, 1), 1);

load gplvmOil

% compute the kernel from results
[K, invK] = computeKernel(X(activeSet, :), theta);

% Visualise the results
gplvmvisualise(X, Y, invK, theta, lbls, meanData, activeSet, 'vectorVisualise', ...
	       'vectorModify');
