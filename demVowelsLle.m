% DEMVOWELSLLE Model the vowels data with a 2-D FGPLVM using RBF kernel.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'vowels';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);
X = lleEmbed(Y, 2);

