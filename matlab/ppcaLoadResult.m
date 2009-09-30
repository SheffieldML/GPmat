function [model, lbls] = ppcaLoadResult(dataSet, number)

% PPCALOADRESULT Load a previously saved result.
% FORMAT
% DESC loads a previously saved PPCA result.
% ARG dataSet : the name of the data set to load.
% ARG number : the number of the PPCA run to load.
% RETURN model : the saved model.
%
% SEEALSO : ppcaLoadResult
%
% COPYRIGHT : Neil D. Lawrence,  2009
  
% MLTOOLS

  [Y, lbls] = lvmLoadData(dataSet);

  dataSet(1) = upper(dataSet(1));
  load(['dem' dataSet 'Ppca' num2str(number)])
  model = ppcaReconstruct(ppcaInfo, Y);
end