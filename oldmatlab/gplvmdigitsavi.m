% GPLVMDIGITSAVI Make AVI files of digits data.

% load data
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

load gplvmDigits1D

[K, invK] = computeKernel(X(activeSet, :), theta);

% Visualise the results
M = gplvmmakeavi1D(X, Y, invK, theta, [], meanData, activeSet, 'imageVisualise', ...
		   'imageModify', 3000, [16 16]);
movie2avi(M, 'digitsFantasy.avi', 'compression', 'none', 'videoname', ...
	  'Fantasy images of Digits', 'FPS', 24)

M = gplvmdatamakeavi1d(X, Y, invK, theta, [], meanData, activeSet, 'imageVisualise', ...
		   'imageModify', [16 16]); 
movie2avi(M, 'digitsData.avi', 'compression', 'none', 'videoname', ...
	  'Data aligned along latent variable', 'FPS', 24)
