function y = gpReversibleDynamicsSamp(model, X);

% GPREVERSIBLEDYNAMICSSAMP Sample from the dynamics for a given input.

% FGPLVM

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
