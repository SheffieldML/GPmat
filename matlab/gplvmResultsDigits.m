% GPLVMRESULTSDIGITS Load and visualise the 2-D results of the digits.

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

load gplvmDigits

[K, invK] = computeKernel(X(activeSet, :), theta);

% Visualise the results
gplvmvisualise(X, Y, invK, theta, lbls, meanData, activeSet, 'imageVisualise', ...
	       'imageModify', [16 16]);
