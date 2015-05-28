function ind = gpBlockIndices(model, blockNo)

% GPBLOCKINDICES Return indices of given block.
% FORMAT
% DESC returns the indices of a given block for the PITC
% approximation.
% ARG model : the model for which the indices are being computed.
% ARG blockNum : the block number for which the indices are
% required.
% RETURN indices : the data indices associated with given block.
%
% SEEALSO : gpComputeAlpha, gpCovGrads, gpLogLikeGradients,
% gpLogLikelihood, gpUpdateAD
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GP


if ~strcmp(model.approx, 'pitc')
  error('Block indices only relevant for pitc approximation.');
else
  if blockNo == 1
    startVal = 1;
  else
    startVal = model.blockEnd(blockNo-1)+1;
  end
  endVal = model.blockEnd(blockNo);
  ind = startVal:endVal;
end
