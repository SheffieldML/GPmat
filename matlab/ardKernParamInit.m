function kern = ardKernParamInit(kern)

% ARDKERNPARAMINIT ARD kernel parameter initialisation.

% KERN


kern.inverseWidth = 1;
kern.rbfVariance = 1;
kern.whiteVariance = 1; 
kern.biasVariance = 1;
kern.linearVariance = 1;
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 5 + kern.inputDimension;


kern.transforms(1).index = [1 2 3 4 5];
kern.transforms(1).type = 'negLogLogit';
kern.transforms(2).index = [6:kern.nParams];
kern.transforms(2).type = 'sigmoid';
