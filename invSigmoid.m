function y=invSigmoid(x)

% INVSIGMOID The inverse of the sigmoid function.
% FORMAT
% DESC returns the inverse of the sigmoid function (which takes the
% form y=log(x/(1-x)). 
% ARG x : the input to the inverse of the sigmoid (should be
% between 0 and 1).
% ARG y : the inverse of the sigmoid.
% 
% SEEALSO : sigmoid
%
% COPYRIGHT : Neil D. Lawrence, 2004

% NDLUTIL

y = log(x./(1-x));
