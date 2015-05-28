function y = sigmoidTransform(x, transform)

% SIGMOIDTRANSFORM Constrains a parameter to be between 0 and 1.
% FORMAT
% DESC contains commands to constrain parameters to be between 0
% and 1 via the sigmoid function.
% ARG x : input argument.
% ARG y : return argument.
% ARG transform : type of transform, 'atox' maps a value into
% the transformed space (i.e. makes it between 0 and 1). 'xtoa' maps the
% parameter back from transformed space to the original
% space. 'gradfact' gives the factor needed to correct gradients
% with respect to the transformed parameter.
% 
% SEEALSO : negLogLogitTransform, expTransform
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2007

% OPTIMI


limVal = 36;
y = zeros(size(x));
switch transform
 case 'atox'
  index = find(x<-limVal);
  y(index) = eps;
  x(index) = NaN;
  index = find(x<limVal);
  y(index) = sigmoid(x(index));
  x(index) = NaN;
  index = find(~isnan(x));
  y(index) = 1-eps;

 case 'xtoa'
  y = invSigmoid(x);
 case 'gradfact'
  y = x.*(1-x);
end
