function noise = noiseReadFromFID(FID)

% NOISEREADFROMFID Load from an FID written by the C++ implementation.
%
%	Description:
%
%	NOISEREADFROMFID(FID, NOISE) reads a noise structure from a file
%	stream created by the C++ implementation of the noises.
%	 Arguments:
%	  FID - the file stream ID.
%	  NOISE - the noise structure created.
%	
%
%	See also
%	MODELREADFROMFID, NOISECREATE, NOISEREADPARAMSFROMFID


%	Copyright (c) 2005, 2008 Neil D. Lawrence


type = readStringFromFID(FID, 'type');
noise = noiseReadParamsFromFID(type, FID);

