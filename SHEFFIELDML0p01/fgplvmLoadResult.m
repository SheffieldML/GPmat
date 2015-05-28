function [model, lbls] = fgplvmLoadResult(dataSet, number)

% FGPLVMLOADRESULT Load a previously saved result.
%
%	Description:
%
%	[MODEL, LBLS] = FGPLVMLOADRESULT(DATASET, NUMBER) loads a previously
%	saved FGPLVM result.
%	 Returns:
%	  MODEL - the saved model.
%	  LBLS - labels of the data set (for visualisation purposes).
%	 Arguments:
%	  DATASET - the name of the data set to load.
%	  NUMBER - the number of the FGPLVM run to load.
%	
%
%	See also
%	FGPLVMLOADRESULT


%	Copyright (c) 2003, 2004, 2005, 2006, 2008 Neil D. Lawrence


[Y, lbls] = lvmLoadData(dataSet);

dataSet(1) = upper(dataSet(1));
load(['dem' dataSet 'Fgplvm' num2str(number)])
model = fgplvmReconstruct(kern, noise, fgplvmInfo, X, Y);