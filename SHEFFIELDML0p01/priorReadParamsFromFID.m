function prior = priorReadParamsFromFID(prior, FID, version)

% PRIORREADPARAMSFROMFID Read prior params from C++ written FID.
%
%	Description:
%
%	PRIOR = PRIORREADPARAMSFROMFID(PRIOR, FID) reads prior parameters
%	from a file written by C++.
%	 Returns:
%	  PRIOR - the prior with the parameters added.
%	 Arguments:
%	  PRIOR - the prior to put the parameters in.
%	  FID - the stream from which parameters are read.
%	
%
%	See also
%	PRIORREADFROMFID, MODELREADFROMFID


%	Copyright (c) 2005, 2006, 2008 Neil D. Lawrence


  if nargin < 3
    version = [];
  end
numParams = readIntFromFID(FID, 'numParams');
params = modelReadFromFID(FID, version);
fhandle = str2func([prior.type 'PriorExpandParam']);
prior = fhandle(prior, params);

