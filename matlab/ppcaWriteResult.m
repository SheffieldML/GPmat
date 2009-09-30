function fileName = ppcaWriteResult(model, dataSet, number)

% PPCAWRITERESULT Write a PPCA result.
% FORMAT
% DESC writes a PPCA result.
% ARG model : the model to write.
% ARG dataSet : the name of the data set to write.
% ARG number : the number of the PPCA run to write.
% RETURN fileName : the file name used to write.
%
% SEEALSO : ppcaLoadResult
%
% COPYRIGHT : Neil D. Lawrence, 2009
  
% PPCA

dataSet(1) = upper(dataSet(1));
type = model.type;
type(1) = upper(type(1));
fileName = ['dem' dataSet type num2str(number)];

ppcaInfo = ppcaDeconstruct(model);

save(fileName, 'ppcaInfo');