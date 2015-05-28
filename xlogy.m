function z = xlogy(x, y)

% XLOGY z = x*log(y) returns zero if x=y=0
% FORMAT 
% DESC computes the function x*log(y) but taking account for the
% fact that the answer is zero if x=y=0.
% ARG x : the x argument in x*log(y).
% ARG y : the y argument in x*log(y).
% RETURN z : returns z where z = x*log(y).
% 
% FORMAT
% DESC computes the function x*log(x), taking account for the fact
% that the answer is zero if x=0.
% ARG x : the argument in x*log(x).
% RETURN y : returns y where y = x*log(x).
%
% COPYRIGHT : Neil D. Lawrence, 2001, 2006

% NDLUTIL

% If there is only one input argument x is assumed to equal y.
if nargin == 1
  y=x;
end
if any(any(x==0))
  if ~issparse(x)
    z = zeros(size(x));
  else
    z = spalloc(size(x, 1), size(x, 2), sum(sum(x~=0)));
  end
  indx = find(x);
  z(indx)= x(indx).*log(y(indx));
else
  z= x.*log(y);
end
