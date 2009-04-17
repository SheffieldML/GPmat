function prior = priorReadFromFID(FID, version)

% PRIORREADFROMFID Read a prior from a C++ written FID.
% FORMAT
% DESC loads in from a file stream the data format produced by 
% C++ implementations.
% ARG FID : the file ID from where the data is loaded.
% RETURN prior : the prior loaded in from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2008
%
% SEEALSO : modelReadFromFID, priorCreate, priorReadParamsFromFID

% PRIOR
if nargin < 2
  version = [];
end
type = readStringFromFID(FID, 'type');
prior = priorCreate(type);
prior = priorReadParamsFromFID(prior, FID, version);


