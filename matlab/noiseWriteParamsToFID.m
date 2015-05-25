function noiseWriteParamsToFID(noise, FID)

% NOISEWRITEPARAMSTOFID Write the noise parameters to a stream.
% FORMAT
% DESC writes noise parameters to  a file stream.
% ARG noise : the noise structure that is being written.
% ARG FID : the file ID of the stream that is being written.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2008
%
% SEEALSO : noiseWriteToFID

% NOISE

writeIntToFID(FID, 'outputDim', noise.numProcess);
writeIntToFID(FID, 'numParams', noise.nParams);
if strcmp(noise.type, 'ncnm') 
  writeIntToFID(FID, 'gammaSplit', noise.gammaSplit);
end
fhandle = str2func([noise.type 'NoiseExtractParam']);
params = fhandle(noise);
doubleMatrixWriteToFID(params, FID);
