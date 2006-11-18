function model = gpDynamicsExpandParam(model, params)

% GPDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP dynamics.
% FORMAT
% DESC takes the given vector of parameters and places them in the
% model structure, it then updates any stored representations that
% are dependent on those parameters, for example kernel matrices
% etc..
% ARG model : the model structure for which parameters are to be
% updated.
% ARG params : a vector of parameters for placing in the model
% structure.
% RETURN model : a returned model structure containing the updated
% parameters.
% 
% SEEALSO : gpExpandParam, gpDynamicsCreate, gpDynamicsExtractParam, modelExtractParam, gpUpdateKernels
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM


% get the current parameter vector
origParam = gpExtractParam(model);

% Substitute for X_u.
origParam(1:model.k*model.q) = params(1:model.k*model.q);

% Subsitute for any parameters to be optimised.
if model.learn
  origParam(model.k*model.q+1:end) = params(model.k*model.q+1:end);
elseif model.learnScales
  switch model.approx
   case 'ftc'
    endVal = length(origParam);
    startVal = endVal-model.d+1;
   case {'dtc', 'fitc', 'pitc'}
    endVal = length(origParam)-1;
    startVal = endVal-model.d+1;
  end
  origParam(startVal:endVal) = params(model.k*model.q+1:end);
end

% Now use the standard gpExpandParam
model = gpExpandParam(model, origParam);