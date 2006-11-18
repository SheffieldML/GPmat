function ind = gpDataIndices(model, dimNo, blockNo)

% GPDATAINDICES Return indices of present data.
% FORMAT
% DESC returns the indices of data which is not missing for a given
% dimension in the GP-LVM.
% ARG model : the model for which the indices are being returned.
% ARG dimNo : the dimension for which the presence of missing data
% is being looked at.
% RETURN ind : indices of training data along that dimension which
% isn't missing.
%
% DESC returns the indices of data which is not missing for a given
% dimension in the GP-LVM and a block number in the PITC approximation.
% ARG model : the model for which the indices are being returned.
% ARG dimNo : the dimension for which the presence of missing data
% is being looked at.
% ARG blockNo : the block number in the PITC approximation for
% which the indices are required.
% RETURN ind : indices of training data along that dimension which
% isn't missing.
%
% SEEALSO : gpCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GP


if nargin > 2
  if ~strcmp(model.approx, 'pitc')
    error('Block number only relevant for pitc approximation.');
  else
    if blockNo == 1
      startVal = 1;
    else
      startVal = model.blockEnd(blockNo-1)+1;
    end
    endVal = model.blockEnd(blockNo);
    if model.isMissingData
      st = min(find(model.indexPresent{dimNo}>=startVal));
      fi = max(find(model.indexPresent{dimNo}<=endVal));
      ind = model.indexPresent{dimNo}(st:fi)';
    else
      ind = startVal:endVal;
    end
  end
else
  if strcmp(model.approx, 'pitc')
    error('Must give block number with PITC approximation');
  else
    if model.isMissingData
      ind = model.indexPresent{dimNo}';
    else
      ind = 1:model.N;
    end
  end
end
