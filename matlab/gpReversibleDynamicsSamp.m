function y = gpReversibleDynamicsSamp(model, X);

% GPREVERSIBLEDYNAMICSSAMP Sample from the dynamics for a given input.
%
% y = gpReversibleDynamicsSamp(model, X);
%

% Copyright (c) 2006 Neil D. Lawrence
% gpReversibleDynamicsSamp.m version 1.1



persistent oldX
if isempty(oldX)
  Xp = [X zeros(size(X))];
else
  Xp = [X X-oldX];
end
oldX = X;
[mu, var] = gpPosteriorMeanVar(model, Xp);
y = gsamp(mu, diag(var), 1);
y = X + y;
