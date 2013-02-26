function z = xlogy(x, y)

% XLOGY z = x*log(y) returns zero if x=y=0
%
%	Description:
%
%	Z = XLOGY(X, Y) computes the function x*log(y) but taking account
%	for the fact that the answer is zero if x=y=0.
%	 Returns:
%	  Z - returns z where z = x*log(y).
%	 Arguments:
%	  X - the x argument in x*log(y).
%	  Y - the y argument in x*log(y).
%
%	Y = XLOGY(X) computes the function x*log(x), taking account for the
%	fact that the answer is zero if x=0.
%	 Returns:
%	  Y - returns y where y = x*log(x).
%	 Arguments:
%	  X - the argument in x*log(x).


%	Copyright (c) 2001, 2006 Neil D. Lawrence


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