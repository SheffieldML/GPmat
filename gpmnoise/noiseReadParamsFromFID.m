function noise = noiseReadParamsFromFID(type, FID)

% NOISEREADPARAMSFROMFID Read the noise parameters from C++ file FID.
% FORMAT
% DESC reads noise parameters from a C++ file stream and returns them in
% the appropriate noise structure.
% ARG type : the type of noise structure that is being read.
% ARG FID : the file ID of the stream that is being read.
% RETURN noise : the noise structure containing the parameters.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2008
%
% SEEALSO : noiseReadFromFID

% NOISE

noise.numProcess = readIntFromFID(FID, 'outputDim');
noise.nParams = readIntFromFID(FID, 'numParams');
noise.type = type;  

if strcmp(type, 'ncnm') 
  noise.gammaSplit = readIntFromFID(FID, 'gammaSplit');
end
noise = noiseParamInit(noise);
params = modelReadFromFID(FID);
fhandle = str2func([noise.type 'NoiseExpandParam']);
noise = fhandle(noise, params);
