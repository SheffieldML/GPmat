% DEMMLPMOUTHDATA2 Try on Ismael's vowels number 2.

lipsScriptLoadBlabla;

Xtemp = acVec;
scales = sqrt(var(mouthPoints));
for i=1:size(mouthPoints, 2)
  mouthPoints(:, i) = mouthPoints(:, i)/scales(i);
end
ytemp = mouthPoints;

x = Xtemp;
t = ytemp;

% Set up network parameters.
nin = size(x, 2);		% Number of inputs.
nhidden = 10;		% Number of hidden units.
nout = size(t, 2);		% Number of outputs.
aw1 = 0.01*ones(1, nin);	% First-layer ARD hyperparameters.
ab1 = 0.01;			% Hyperparameter for hidden unit biases.
aw2 = 0.01;			% Hyperparameter for second-layer weights.
ab2 = 0.01;			% Hyperparameter for output unit biases.
beta = 50.0;			% Coefficient of data error.

% Create and initialize network.
prior = mlpprior(nin, nhidden, nout, aw1, ab1, aw2, ab2);
% Create and initialize network weight vector.
net = mlp(nin, nhidden, nout, 'linear', prior, beta);

% Set up vector of options for the optimiser.
nouter = 3;			% Number of outer loops.
ninner = 10;			% Number of innter loops.
options = zeros(1,18);		% Default options vector.
options(1) = 1;			% This provides display of error values.
options(2) = 1.0e-7;		% Absolute precision for weights.
options(3) = 1.0e-7;		% Precision for objective function.
options(14) = 1000;		% Number of training cycles in inner loop. 

% Train using scaled conjugate gradients, re-estimating alpha and beta.
for k = 1:nouter
  net = netopt(net, options, x, t, 'scg');
  [net, gamma] = evidence(net, x, t, ninner);
  fprintf(1, '\nRe-estimation cycle %d:\n', k);
  fprintf(1, '  alpha =  %8.5f\n', net.alpha);
  fprintf(1, '  beta  =  %8.5f\n', net.beta);
  fprintf(1, '  gamma =  %8.5f\n\n', gamma);
end

fprintf(1, 'true beta: %f\n', 1/(noise*noise));
