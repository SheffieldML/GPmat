function model = expvarMeanExpandParam(model, params)

% EXPVARMEANEXPANDPARAM Update expvarMean mapping with new vector of parameters.
% FORMAT
% DESC takes a vector of EXPVARMEAN parameters and places them in their
% respective positions in the EXPVARMEAN model. The function is a wrapper
% for the mlpunpa command.
% ARG model : the model in which the weights are to be placed.
% ARG params : a vector of the weights to be placed in the model.
% RETURN model : the model with the weights distributed in the
% correct places.
%
% SEEALSO : expvarMeanCreate, expvarMeanExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

% simply update the parameters of the kernel.
model.kern = kernExpandParam(model.kern, params);
