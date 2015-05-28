function y = negLogLogit(x)

% NEGLOGLOGIT Function which returns the negative log of the logistic function.
% FORMAT
% DESC computes the negative log of the logistic (sigmoid)
% function, which is also the integral of the sigmoid function.
% ARG x : input locations.
% RETURN y : the negative log of the logistic.
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% SEEALSO : sigmoid

% NDLUTIL

y = log(1+exp(x));
