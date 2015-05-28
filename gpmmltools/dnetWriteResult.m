function fileName = dnetWriteResult(model, dataSet, number)

% DNETWRITERESULT Write a DNET result.
% FORMAT
% DESC writes a DNET result.
% ARG model : the model to write.
% ARG dataSet : the name of the data set to write.
% ARG number : the number of the DNET run to write.
% RETURN fileName : the file name used to write.
%
% SEEALSO : dnetLoadResult
%
% COPYRIGHT : Neil D. Lawrence, 2009
  
% MLTOOLS

dataSet(1) = upper(dataSet(1));
type = model.type;
type(1) = upper(type(1));
fileName = ['dem' dataSet type num2str(number)];

[mapping, dnetInfo] = dnetDeconstruct(model);

save(fileName, 'mapping', 'dnetInfo');
