function priorWriteToFID(FID, prior)

% PRIORWRITETOFID Write a prior to a C++ stream.
% FORMAT
% DESC writes to a file stream a prior.
% ARG FID : the file ID from where the data is loaded.
% ARG prior : the prior loaded in from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : modelWriteToFID, priorCreate, priorWriteParamsToFID

% PRIOR

writeVersionToFID(FID, 0.2);
writeStringToFID(FID, 'baseType', 'prior');
writeStringToFID(FID, 'type', prior.type);
priorWriteParamsToFID(prior, FID);


