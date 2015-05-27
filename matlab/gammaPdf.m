function [p, lnp] = gammaPdf(x, a, b)

% GAMMAPDF PDF for the Gamma distribution.
% FORMAT
% DESC computes the pdf of the gamma distribution.
% ARG x : locations where the PDF is to be computed.
% ARG a : shape parameter of the gamma distribution.
% ARG b : rate parameter of the gamma distribuion (inverse scale).
% RETURN p : probability of the gamma distribution.
% RETURN lnp : log of the probability.
%
% SEEALSO cumGamma
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% NDLUTIL


lnp = a*log(b) - gammaln(a) + (a-1)*log(x) - b*x;
p = exp(lnp);
