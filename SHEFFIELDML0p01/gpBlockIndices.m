function ind = gpBlockIndices(model, blockNo)

% GPBLOCKINDICES Return indices of given block.
%
%	Description:
%
%	INDICES = GPBLOCKINDICES(MODEL, BLOCKNUM) returns the indices of a
%	given block for the PITC approximation.
%	 Returns:
%	  INDICES - the data indices associated with given block.
%	 Arguments:
%	  MODEL - the model for which the indices are being computed.
%	  BLOCKNUM - the block number for which the indices are required.
%	gpLogLikelihood, gpUpdateAD
%	
%
%	See also
%	GPCOMPUTEALPHA, GPCOVGRADS, GPLOGLIKEGRADIENTS, 


%	Copyright (c) 2006 Neil D. Lawrence



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
