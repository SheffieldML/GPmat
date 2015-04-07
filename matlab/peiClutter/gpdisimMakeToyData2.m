function [genes, params, proteins] = gpdisimMakeToyData2(numGenes, numPoints, numReplicates)

times = (1:numPoints)';

delta = 1;
sigma = gamrnd(1, 10);
D = gamrnd(ones(1, numGenes), 10);
S = gamrnd(ones(1, numGenes), 10);
B = gamrnd(ones(1, numGenes), 10);
l = 3;


t0 = (1:.01:numPoints)';
k_in = kernCreate(t0, 'rbf');
k_in = kernExpandParam(k_in, [log(2/l^2), 0]);
k = kernCompute(k_in, t0);
y0 = log(1 + exp(gsamp(zeros(size(t0)), k, 1)))';

for k = 1:numReplicates,

  y = y0 + .1 * randn(size(y0));
  [t_f, f] = ode45(@(t, x) sigma * interp1(t0, y, t) - delta * x,...
		   [1, numPoints], 0);

  z = interp1(t0, y, times);
  proteins{k} = interp1(t_f, f, times);
  
  for j = 1:numGenes,
    [t_x, x_j] = ode45(@(t, x) B(j) + S(j) * interp1(t_f, f, t) - D(j) * x,...
		       [1, numPoints], 0);
    
    z = [z, interp1(t_x, x_j, times)];
  end

  genes{k} = z + .1 * randn(size(z));
  %genes{k} = reshape(genes{k}, [numPoints, numGenes+1]);
end

params = struct('sigma', sigma, 'delta', delta, 'D', D, 'S', S, 'B', B, 'l', l);
