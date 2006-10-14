function ind = gpDataIndices(model, dimNo, blockNo)

% GPDATAINDICES Return indices of present data.

% FGPLVM


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
