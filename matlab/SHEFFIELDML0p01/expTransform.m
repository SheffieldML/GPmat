function y = expTransform(x, transform)

% EXPTRANSFORM Constrains a parameter to be positive through exponentiation.
%
%	Description:
%
%	EXPTRANSFORM(X, Y, TRANSFORM) contains commands to constrain
%	parameters to be positive via exponentiation.
%	 Arguments:
%	  X - input argument.
%	  Y - return argument.
%	  TRANSFORM - type of transform, 'atox' maps a value into the
%	   transformed space (i.e. makes it positive). 'xtoa' maps the
%	   parameter back from transformed space to the original space.
%	   'gradfact' gives the factor needed to correct gradients with
%	   respect to the transformed parameter.
%	
%
%	See also
%	NEGLOGLOGITTRANSFORM, SIGMOIDTRANSFORM


%	Copyright (c) 2004, 2005, 2006, 2007 Neil D. Lawrence



limVal = 36;
y = zeros(size(x));
switch transform
 case 'atox'
  index = find(x<-limVal);
  y(index) = exp(-limVal);
  x(index) = NaN;
  index = find(x<limVal);
  y(index) = exp(x(index));
  x(index) = NaN;
  index = find(~isnan(x));
  if ~isempty(index)
    y(index) = exp(limVal);
  end
 case 'xtoa'
  y = log(x);
 case 'gradfact'
  y = x;
end
