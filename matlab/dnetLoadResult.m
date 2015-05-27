function [model, lbls] = dnetLoadResult(dataSet, number)

% DNETLOADRESULT Load a previously saved result.
% FORMAT
% DESC loads a previously saved DNET result.
% ARG dataSet : the name of the data set to load.
% ARG number : the number of the DNET run to load.
% RETURN model : the saved model.
%
% SEEALSO : dnetLoadResult
%
% COPYRIGHT : Neil D. Lawrence,  2009
  
% MLTOOLS

  [Y, lbls] = lvmLoadData(dataSet);

  dataSet(1) = upper(dataSet(1));
  load(['dem' dataSet 'Dnet' num2str(number)])
  model = dnetReconstruct(mapping, dnetInfo, Y);
end
