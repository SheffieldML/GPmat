% GPLVMRESULTSBRENDAN Load and visualise the 2-D results of Brendan's face.

% load data
load frey_rawface.mat

Y = double(ff)';

meanData = zeros(1, size(Y, 2));

load gplvmBrendan

[K, invK] = computeKernel(X(activeSet, :), theta);

% Visualise the results
gplvmvisualise(X, Y, invK, theta, [], meanData, activeSet, 'imageVisualise', ...
	       'imageModify', [20 28]);
