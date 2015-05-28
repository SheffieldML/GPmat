function y = negLogLogitTransform(x, transform)

% NEGLOGLOGITTRANSFORM Constrains a parameter to be positive.
% FORMAT
% DESC contains commands to constrain parameters to be positive via
% log(1+exp(x)).
% ARG x : input argument.
% ARG y : return argument.
% ARG transform : type of transform, 'atox' maps a value into
% the transformed space (i.e. makes it positive). 'xtoa' maps the
% parameter back from transformed space to the original
% space. 'gradfact' gives the factor needed to correct gradients
% with respect to the transformed parameter.
% 
% SEEALSO : expTransform, sigmoidTransform
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2007

% OPTIMI


y = zeros(size(x));
limVal = 36;
switch transform
 case 'atox'
  index = find(x<-limVal);
  y(index) = eps;
  x(index) = NaN;
  index = find(x<limVal);
  y(index) = log(1+exp(x(index)));
  x(index) = NaN;
  index = find(~isnan(x));
  y(index) = x(index);
 case 'xtoa'
  index = find(x<limVal);
  y(index) = log(exp(x(index))-1);
  index = find(x>=limVal);
  y(index) = x(index);
 case 'gradfact'
  index = find(x>limVal);
  y(index) = 1;
  index = find(x<=limVal);
  y(index) = 1-exp(-x(index));
end
  
