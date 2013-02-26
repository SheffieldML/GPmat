function prior = priorReadFromFID(FID, version)

% PRIORREADFROMFID Read a prior from a C++ written FID.
%
%	Description:
%
%	PRIOR = PRIORREADFROMFID(FID) loads in from a file stream the data
%	format produced by C++ implementations.
%	 Returns:
%	  PRIOR - the prior loaded in from the file.
%	 Arguments:
%	  FID - the file ID from where the data is loaded.
%	
%
%	See also
%	MODELREADFROMFID, PRIORCREATE, PRIORREADPARAMSFROMFID


%	Copyright (c) 2005, 2006, 2008 Neil D. Lawrence

if nargin < 2
  version = [];
end
type = readStringFromFID(FID, 'type');
prior = priorCreate(type);
prior = priorReadParamsFromFID(prior, FID, version);


