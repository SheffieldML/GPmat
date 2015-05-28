function ll = modelPointLogLikelihood(model, varargin)

% MODELPOINTLOGLIKELIHOOD Compute the log likelihood of a given point.
% FORMAT
% DESC computes the log likelihood of the given model.
% ARG model : the model for which the log likelihood is to be
% computed.
% ARG P1, P2 ... : additional arguments as required.
% RETURN ll : the log likelihood of the given data point.
%
% SEEALSO : modelLogLikelihood, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

fhandle = str2func([model.type 'PointLogLikelihood']);
ll = fhandle(model, varargin{:});
