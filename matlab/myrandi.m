function R = myrandi(from, to, N, D)

% MYRANDI Create N D-dimensional random numbers from 'from' to 'to'.
% DESC: Should be similar to MATLAB function 'randi' which is not present in older distributions.
%
% COPYRIGHT: Andreas C. Damianou, 2013
% SHEFFIELDML

if nargin < 4, D = 1; end
if nargin < 3,  N = 1; end
if nargin < 2, error('At least two arguments needed'); end

R = zeros(N,D);

for d = 1:D
    R(:,d) = from + (rand(N,1) * (to-from));
end
