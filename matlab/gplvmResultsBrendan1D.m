% GPLVMRESULTSBRENDAN1D Load and visualise the 1-D results for the faces.

% load data
load frey_rawface.mat

Y = double(ff)';

meanData = mean(Y);
Y = Y  - repmat(meanData, size(Y, 1), 1);

load gplvmBrendan1D.mat

[K, invK] = computeKernel(X(activeSet, :), theta);

% Visualise the results
gplvmvisualise1D(X, Y, invK, theta, [], meanData, activeSet, 'imageVisualise', ...
		 'imageModify', [20 28]);

