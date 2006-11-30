function model = gpsimMapUpdateG(model)

% GPSIMMAPUDATEG Update the nonlinear transformation of f.
% FORMAT
% DESC updates the fields of model associated with the non-linear
% transformation of f, namely g, g_grad and g_grad2.
% ARG model : the model to be updated.
% ARG model : the model with the updated g representation.
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% SEEALSO : gpsimMapFunctionalExpandParam, gpsimMapCreate

% GPSIM

switch model.nonLinearity
 case 'linear'
  model.g=model.f;        %linear case
  model.g_grad=ones(size(model.g));
  model.g_grad2=zeros(size(model.g));
 case 'exp'
  model.g=exp(model.f);   %positive TF concentrations
  model.g_grad=model.g;
  model.g_grad2=model.g;
 otherwise
  error('Invalid non-linearity.')
end

