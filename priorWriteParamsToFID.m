function prior = priorWriteParamsToFID(prior, FID)

% PRIORWRITEPARAMSTOFID Write prior params from C++ written FID.
% FORMAT
% DESC writes prior parameters from a file written by C++.
% ARG prior : the prior to put the parameters in.
% ARG FID : the stream from which parameters are write.
% RETURN prior : the prior with the parameters added.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : priorWriteToFID, modelWriteToFID

% PRIOR

writeIntToFID(FID, 'numParams', prior.nParams);
fhandle = str2func([prior.type 'PriorExtractParam']);
params = fhandle(prior);
doubleMatrixWriteToFID(params, FID);

