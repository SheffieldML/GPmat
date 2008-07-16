function model = gpReversibleDynamicsCreate(q, d, latentVals, options)

% GPREVERSIBLEDYNAMICSCREATE Create the dynamics model. 

% FGPLVM

if nargin < 4
  options = gpReversibleDynamicsOptions('ftc');
end

diffs = latentVals(2:end, :) - latentVals(1:end-1, :);
X = [latentVals(2:end-1, :) diffs(1:end-1, :)];
y = diffs(2:end, :);
model = gpCreate(2*q, d, X, y, options);
model.learn = 0;
model.type = 'gpReversibleDynamics';

