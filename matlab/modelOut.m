function [Y, Phi] = modelOut(model, X, varargin)

% MODELOUT Give the output of a model for given X.
% FORMAT
% DESC gives the output of the model for a given input X. For
% latent variable models it gives a position in data space given a
% position in latent space.
% ARG model : structure specifying the model.
% ARG X : input location(s) for which output is to be computed.
% RETURN Y : output location(s) corresponding to given input
% locations.
%
% FORMAT
% DESC gives the output of the model for a given input X. For
% latent variable models it gives a position in data space given a
% position in latent space.
% ARG model : structure specifying the model.
% ARG X : input location(s) for which output is to be computed.
% RETURN Phi : output basis function(s) corresponding to given input
% RETURN Y : output location(s) corresponding to given input
% locations.
%
% SEEALSO : modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Cark Henrik Ek, 2008

% MLTOOLS

fhandle = str2func([model.type 'Out']);
if nargout > 1
  [Y, Phi] = fhandle(model, X, varargin{:});
else
  Y = fhandle(model, X, varargin{:});
end
if(isfield(model,'indexOut')&&~isempty(model.indexOut))
  Y(:,setdiff(1:1:size(Y,2),model.indexOut)) = NaN;
end
