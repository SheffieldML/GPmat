function componentKernWriteParamsToFID(kern, FID)

% COMPONENTKERNWRITEPARAMSTOFID Write a component based kernel to a stream.
%
%	Description:
%
%	COMPONENTKERNWRITEPARAMSTOFID(KERN, FID) writes the components from
%	a component kernel to a stream.
%	 Arguments:
%	  KERN - the base kernel the components come from.
%	  FID - the output file stream to write to.
%	
%
%	See also
%	MODELWRITETOFID, KERNWRITETOFID


%	Copyright (c) 2008 Neil D. Lawrence


writeIntToFID(FID, 'inputDim', kern.inputDimension);
writeIntToFID(FID, 'numParams', kern.nParams);
writeIntToFID(FID, 'numKerns', length(kern.comp));

for i=1:length(kern.comp)
  kernWriteToFID(kern.comp{i}, FID);
end

