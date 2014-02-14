function y = doTransform(x, type, transform)

% DOTRANSFORM Wrapper to apply a transform
%
% FORMAT
%
% DESC Wrapper to apply a transform, possibly with futher parameters
%
% ARG x : input argument.
%
% ARG type : type of transform, 'atox' maps a value into
% the transformed space (i.e. makes it between A and B). 'xtoa' 
% maps the parameter back from transformed space to the original
% space. 'gradfact' gives the factor needed to correct gradients
% with respect to the transformed parameter, that is, it gives
% the gradient dx/da where x is the transformed parameter and a
% is the untransformed parameter.
%
% ARG transform: a structure with fields 'type' to indicate the
% type of transfomation (e.g. 'exp', 'identity', 'sigmoid' or
% 'sigmoidab') and optionally 'transformsettings' for further
% parameters
%
% OUTPUT y : return argument.
% 
% SEEALSO : negLogLogitTransform, expTransform
%
% COPYRIGHT : Antti Honkela, 2013

% SHEFFIELDML

fhandle = str2func([transform.type 'Transform']);
if isfield(transform, 'transformsettings'),
  y = fhandle(x, type, transform.transformsettings);
else
  y = fhandle(x, type);
end
