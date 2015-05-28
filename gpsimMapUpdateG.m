function model = gpsimMapUpdateG(model)

% GPSIMMAPUPDATEG Update the nonlinear transformation of f.
% FORMAT
% DESC updates the fields of model associated with the non-linear
% transformation of f, namely g, g_grad and g_grad2.
% ARG model : the model to be updated.
% ARG model : the model with the updated g representation.
%
% COPYRIGHT : Neil D. Lawrence, 2006
%  
% MODIFIED : Pei Gao, 2008
%
% SEEALSO : gpsimMapFunctionalExpandParam, gpsimMapCreate

% SHEFFIELDML

% Remove mean function value from m (if mean function present).
if isfield(model, 'meanFunction') & ~isempty(model.meanFunction)
  f = model.f + modelOut(model.meanFunction, model.mapt);
else
  f = model.f;
end

% make sure f is a column vector
f = reshape(f, length(f), 1);

if model.ngParam == 0

switch model.nonLinearity
   case 'linear'
    model.g=f;        %linear case  
    model.g_grad=ones(size(model.g));
    model.g_grad2=zeros(size(model.g));
    model.g_grad3=zeros(size(model.g));
    model.isConcave = true; % Is the log likelihood concave?
   case 'exp'
    model.g=exp(f);   %positive TF concentrations
%   model.gParam
%   if ~isempty(model.gParam)
%   model.ngParam = length(model.gParam(1,:));
%   model.dg = [];
%   model.dg2 = [];
%   end  
    model.g_grad=model.g;
    model.g_grad2=model.g;
    model.g_grad3=model.g;
    model.isConcave = true;
   case 'quadratic'
    model.g=f.*f;
    model.g_grad=2*f;
    model.g_grad2=2*ones(size(f));
    model.isConcave = false;
   case 'negLogLogit'
    model.g=log(1+exp(f));   %positive TF concentrations
    model.g_grad=sigmoid(f);
    model.g_grad2=model.g_grad.*(1-model.g_grad);
    model.isConcave = true;
   case 'sigmoid'
    model.g = sigmoid(f);
    model.g_grad = model.g.*(1 - model.g);
    model.g_grad2 = model.g.*(1-model.g) - 2*model.g.*(model.g.*(1-model.g));
    model.isConcave = false;
   case 'repression'
    model.g = gpsimModelFunctions(model,'g');
    model.g_grad = gpsimModelFunctions(model,'grad');
    model.g_grad2 = gpsimModelFunctions(model,'grad2');
    model.g_grad3 = gpsimModelFunctions(model,'grad3');;
    model.isConcave = true; %% ??
   case 'activation'
    model.g = gpsimModelFunctions(model,'g');
    model.g_grad = gpsimModelFunctions(model,'grad');
    model.g_grad2 = gpsimModelFunctions(model,'grad2');
    model.g_grad3 = gpsimModelFunctions(model,'grad3');
    model.isConcave = true; %% ??
    if ~isempty(model.gParam)
      model.ngParam = length(model.gParam(1,:));
      model.dg = gpsimModelFunctions(model,'paramGrad');
      model.dg2 = gpsimModelFunctions(model,'paramGrad2');
    end  
   otherwise
    error('Invalid non-linearity.')
end

else
  model.g = gpsimModelFunctions(model,'g');
  model.g_grad = gpsimModelFunctions(model,'grad');
  model.g_grad2 = gpsimModelFunctions(model,'grad2');
  model.g_grad3 = gpsimModelFunctions(model,'grad3');
  model.isConcave = true;
  model.dg = gpsimModelFunctions(model,'paramGrad');
  model.dggrad = gpsimModelFunctions(model,'paramGgrad');
  model.dggrad2 = gpsimModelFunctions(model,'paramGgrad2');
end


