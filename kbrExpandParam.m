function model = kbrExpandParam(model, params,dim);

% KBREXPANDPARAM Create model structure from KBR model's parameters.
% FORMAT
% DESC takes a vector of KBR weights and centres and places them in their
% respective positions in the KBR model. 
% ARG model : the model in which the weights are to be placed.
% ARG params : a vector of the weights to be placed in the model.
% RETURN model : the model with the weights distributed in the
% correct places.
%
% SEEALSO : kbrunpak, kbrCreate, kbrExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% MODIFICATIONS : Carl Henrik Ek, 2008
%
% MLTOOLS

if(nargin<3)
  startVal = 1;
  endVal = model.numData*model.outputDim;
  model.A = reshape(params(1:endVal), model.numData, model.outputDim);
  model.bias = params(endVal+1:end);
else
  model.A(:,dim) = params(1:1:model.numData*length(dim));
  model.bias(dim) = params(model.numData*length(dim)+1:1:end);
end
