function noiseWriteParamsToFID(noise, FID)

% NOISEWRITEPARAMSTOFID Write the noise parameters to a stream.
%
%	Description:
%
%	NOISEWRITEPARAMSTOFID(NOISE, FID) writes noise parameters to  a file
%	stream.
%	 Arguments:
%	  NOISE - the noise structure that is being written.
%	  FID - the file ID of the stream that is being written.
%	
%
%	See also
%	NOISEWRITETOFID


%	Copyright (c) 2005, 2008 Neil D. Lawrence


writeIntToFID(FID, 'outputDim', noise.numProcess);
writeIntToFID(FID, 'numParams', noise.nParams);
if strcmp(noise.type, 'ncnm') 
  writeIntToFID(FID, 'gammaSplit', noise.gammaSplit);
end
fhandle = str2func([noise.type 'NoiseExtractParam']);
params = fhandle(noise);
doubleMatrixWriteToFID(params, FID);
