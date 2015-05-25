% GPLVMBRENDANAVI Make AVI files of faces data.

% load data
load frey_rawface.mat

Y = double(ff)';

meanData = mean(Y);
Y = Y  - repmat(meanData, size(Y, 1), 1);

load gplvmBrendan1D.mat

[K, invK] = computeKernel(X(activeSet, :), theta);

% Visualise the results
M = gplvmmakeavi1D(X, Y, invK, theta, [], meanData, activeSet, 'imageVisualise', ...
		   'imageModify', 1965, [20 28]);
movie2avi(M, 'brendanFantasy.avi', 'compression', 'none', 'videoname', ...
	  'Fantasy images of Face', 'FPS', 24)

M = gplvmdatamakeavi1d(X, Y, invK, theta, [], meanData, activeSet, 'imageVisualise', ...
		   'imageModify', [20 28]); 
movie2avi(M, 'brendanData.avi', 'compression', 'none', 'videoname', ...
	  'Data aligned along latent variable', 'FPS', 24)
