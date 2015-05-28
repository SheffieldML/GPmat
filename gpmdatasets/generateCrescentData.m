function [X, y] = generateCrescentData(numDataPart);

% GENERATECRESCENTDATA Generate crescent data.
% FORMAT
% DESC generates data in a crescent shape.
% ARG numDataPart : number of data in each part of the shape
% (there are four parts, so you get 4 times this amount of data).
% RETURN X : input data locations.
% RETURN y : target data locations.
%
% SEEALSO : mapLoadData
%
% COPYRIGHT : Neil D. Lawrence, 2004

% DATASETS

% Generate a toy data-set
R = [sqrt(2)/2 -sqrt(2)/2; sqrt(2)/2 sqrt(2)/2];
D = [2 0; 0 1];
meanOne = [4 4];
meanTwo = [0 4];
meanThree = -[4 4];
meanFour = -[0 4];

X1 = [randn(numDataPart,2)*D*R]  - repmat(meanOne, numDataPart, 1);
X2 = [randn(numDataPart,2)*D*R'] - repmat(meanTwo, numDataPart, 1);
X3 = [randn(numDataPart,2)*D*R]  - repmat(meanThree, numDataPart, 1);
X4 = [randn(numDataPart,2)*D*R'] - repmat(meanFour, numDataPart, 1);

X = [X1; X2; X3; X4];
y = [ones(2*numDataPart, 1); -ones(2*numDataPart, 1)];
