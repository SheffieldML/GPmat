function componentKernWriteParamsToFID(kern, FID)

% COMPONENTKERNWRITEPARAMSTOFID Write a component based kernel to a stream.
% FORMAT
% DESC writes the components from a component kernel to a stream.
% ARG kern : the base kernel the components come from.
% ARG FID : the output file stream to write to.
%
% SEEALSO : modelWriteToFID, kernWriteToFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% KERN

writeIntToFID(FID, 'inputDim', kern.inputDimension);
writeIntToFID(FID, 'numParams', kern.nParams);
writeIntToFID(FID, 'numKerns', length(kern.comp));

for i=1:length(kern.comp)
  kernWriteToFID(kern.comp{i}, FID);
end

