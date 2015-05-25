function [outputs, params] = gpdisimSimulateSystem(input, params, times),

N = size(input, 2);

if ~isstruct(params),
  numGenes = params;
  params = [];
  params.delta = 1;
  params.sigma = gamrnd(1, 10);
  params.D = gamrnd(ones(1, numGenes), 10);
  params.S = gamrnd(ones(1, numGenes), 10);
  params.B = gamrnd(ones(1, numGenes), 10);
end

if nargin < 3,
  times = (1:12)';
end

numGenes = length(params.D);

sigma = params.sigma;
delta = params.delta;
B = params.B;
D = params.D;
S = params.S;

for k = 1:N,
  y = input(:, k);
  [t_f, f] = ode45(@(t, x) sigma * interp1(times, y, t) - delta * x,...
		   [min(times), max(times)], 0);

  z = y;
  proteins{k} = interp1(t_f, f, times);
  
  for j = 1:numGenes,
    [t_x, x_j] = ode45(@(t, x) B(j) + S(j) * interp1(t_f, f, t) - D(j) * x,...
		       [min(times), max(times)], 0);
    
    z = [z, interp1(t_x, x_j, times)'];
  end

  outputs{k} = z;
  %genes{k} = reshape(genes{k}, [numPoints, numGenes+1]);
end
