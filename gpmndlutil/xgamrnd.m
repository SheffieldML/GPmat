function r = xgamrnd(a, b),
% XGAMRND Draw a sample from the gamma distribution.
% FORMAT
% DESC Generates a random draw from the gamma distribution.
%   The function uses both Matlab's RAND and RANDN generators, so to
%   seed it you must seed both RAND and RANDN.
% ARG a : shape parameter of the gamma distribution.
% ARG b : rate parameter of the gamma distribuion (inverse scale).
% RETURN r : a draw of random numbers from the gamma distribution.
%   If a and b are matrices, r will have the same shape.
%   If either is scalar, r will have the same shape as the one that
%   is a matrix.
%
% SEEALSO rand, randn
%
% COPYRIGHT : Antti Honkela, 2010
  
% NDLUTIL

% Based on Marsaglia and Tsang, "A Simple Method for generating gamma
% variables", ACM Transactions on Mathematical Software 26(3):363-372
% (2000).

% Implemented in Mex.
