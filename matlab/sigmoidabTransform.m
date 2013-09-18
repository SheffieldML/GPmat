function y = sigmoidabTransform(x, transform, transformsettings)

% SIGMOIDABTRANSFORM Constrains a parameter to be between A and B
% by a scaled logistic sigmoid function.
%
% FORMAT
%
% DESC contains commands to constrain parameters to be between A
% and B via the sigmoid function, y=A+(B-A)/(1+exp(-x)).
%
% ARG x : input argument.
%
% ARG transform : type of transform, 'atox' maps a value into
% the transformed space (i.e. makes it between A and B). 'xtoa' 
% maps the parameter back from transformed space to the original
% space. 'gradfact' gives the factor needed to correct gradients
% with respect to the transformed parameter, that is, it gives
% the gradient dx/da where x is the transformed parameter and a
% is the untransformed parameter.
%
% ARG transformsettings (first element of varargin): vector [A B] 
% giving the minimum and maximum values A and B for the transformed 
% parameter. If not given, assume A=0 and B=1 as in the function 
% sigmoidTransform.
%
% OUTPUT y : return argument.
% 
% SEEALSO : negLogLogitTransform, expTransform
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2007
%
% COPYRIGHT : Jaakko Peltonen, 2011
%
% COPYRIGHT : Antti Honkela, 2012

% OPTIMI


% compute logistic sigmoid between A and B
  
A=transformsettings(1);
B=transformsettings(2);

limVal = 36;

y = zeros(size(x));
switch transform
  case 'atox'
    I1 = x < -limVal;
    y(I1) = A+(B-A)*eps;

    I2 = x > limVal;
    y(I2) = A+(B-A)*(1-eps);
    
    I3 = ~I1 & ~I2;
    y(I3) = A+(B-A)*sigmoid(x(I3));
    
  case 'xtoa'
    minval_sigmoid=A+(B-A)*eps;
    maxval_sigmoid=A+(B-A)*(1-eps);

    I1 = x<=minval_sigmoid;
    y(I1) = -limVal;
    
    I2 = x>=maxval_sigmoid;
    y(I2) = limVal;

    I3 = ~I1 & ~I2;
    y(I3) = invSigmoid((x(I3)-A)/(B-A));
    
  case 'gradfact'
    y = (x-A).*(B-x)/(B-A);
end;
