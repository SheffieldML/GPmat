function fileName = fgplvmWriteResult(model, dataSet, number)

% FGPLVMWRITERESULT Write a FGPLVM result.
%
%	Description:
%
%	FILENAME = FGPLVMWRITERESULT(MODEL, DATASET, NUMBER) writess an
%	FGPLVM result.
%	 Returns:
%	  FILENAME - the file name used to write.
%	 Arguments:
%	  MODEL - the model to write.
%	  DATASET - the name of the data set to write.
%	  NUMBER - the number of the FGPLVM run to write.
%	
%
%	See also
%	FGPLVMLOADRESULT


%	Copyright (c) 2009, 2011 Neil D. Lawrence


dataSet(1) = upper(dataSet(1));
type = model.type;
type(1) = upper(type(1));
fileName = ['dem' dataSet type num2str(number)];

[kern, noise, fgplvmInfo, X] = fgplvmDeconstruct(model);

save(fileName, 'kern', 'noise', 'fgplvmInfo', 'X');