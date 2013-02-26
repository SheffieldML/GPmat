function fileName = gpWriteResult(model, dataSet, number)

% GPWRITERESULT Write a GP result.
%
%	Description:
%
%	FILENAME = GPWRITERESULT(MODEL, DATASET, NUMBER) writes a GP result.
%	 Returns:
%	  FILENAME - the file name used to write.
%	 Arguments:
%	  MODEL - the model to write.
%	  DATASET - the name of the data set to write.
%	  NUMBER - the number of the GP run to write.
%	
%
%	See also
%	GPLOADRESULT


%	Copyright (c) 2009 Neil D. Lawrence


dataSet(1) = upper(dataSet(1));
type = model.type;
type(1) = upper(type(1));
fileName = ['dem' dataSet type num2str(number)];

[kern, noise, gpInfo] = gpDeconstruct(model);

save(fileName, 'kern', 'noise', 'gpInfo');