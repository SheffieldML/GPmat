function [r, lls, fitparams] = drosMLRank(drosexp, drosTF, tf, probes, t, do_plots),

% DROSODERANK Rank targets by fitting an ODE to raw expression data.
% FORMAT
% DESC Rank targets by fitting an ODE to raw expression data.
% ARG drosexp : the drosexp data structure from drosLoadData
% ARG drosTF : the drosTF data structure from drosLoadData
% ARG tf : TF symbol, should be in {'bap', 'bin', 'mef2', 'tin', 'twi'}
% ARG probes : probe sets to rank
% ARG t : time points to use (optional, default=all)
% RETURN ranking : indeces of genes in drosexp in ranked order
%
% SEEALSO : drosLoadData, drosScoreTFTargetList
%
% COPYRIGHT : Antti Honkela, 2010

% SHEFFIELDML

if nargin < 5,
  t = 1:12;
end
if nargin < 6,
  do_plots = 0;
end
N = length(drosexp.probes);

[tfprof, tfvar, gene, times, scale, rawExp, rawVar] = drosGetData(drosexp, {drosTF.probes.(tf)});

tfprof = cat(2, tfprof{:});
tfvar = cat(2, tfvar{:});
tfprof = tfprof(t, :);
tfvar = tfvar(t, :);

%params = [0, 0, 0, 0];
%m = integrateODE(tfprof, t, params);

options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);
%options2 = optimset('MaxFunEvals', 100000, 'MaxIter', 100000);
options2=options;

fitparams = {};
lls = zeros(size(probes));

for k=1:length(probes),
  fprintf('Scoring %s (%d/%d)...\n', probes{k}, k, length(probes));
  [y, yvar, gene, times, scale, rawExp, rawVar] = drosGetData(drosexp, probes(k));
  y = cat(2, y{:});
  yvar = cat(2, yvar{:});
  y = y(t, :);
  yvar = yvar(t, :);

  [fitparams{k}, lls(k)] = fminsearch(@(params) ODEObjective(tfprof, tfvar, y, yvar, t, mean(tfprof, 2), params), [0, 0, 0, -2], options);
  %[fitparams{k}, lls(k)] = fminsearch(@(params) ODEObjective(tfprof, tfvar, y, yvar, t, params(1:12), params(13:16)), [mean(tfprof, 2)', fitparams{k}], options2);
  if do_plots,
    plot(t, tfprof, 'b'); hold on
    plot(t, y, 'r');
    plot(t, integrateODE(mean(tfprof, 2), t, fitparams{k}), 'g'); hold off;
    %plot(t, integrateODE(fitparams{k}(1:12)', t, fitparams{k}(13:16)), 'g'); hold off;
    pause
  end
end

lls = -lls;
I = drosFindGeneinds(drosexp, probes, 0, 1);
[foo, J] = sort(lls, 'descend');
r = drosRemoveDuplicateGenes(drosexp, I(J));
%r = gaussianLogLikelihood(m, y, yvar);


function m = integrateODE(f, t, params),

Nreps = size(f, 2);
t_ext = [0, t];
t_ext = t_ext(ones(1, Nreps), :)';

f_ext = [zeros(1, Nreps); f];

mu = exp(params(1));
D = exp(params(2));
S = exp(params(3));
delta = exp(params(4));

C = zeros(length(t)+1, Nreps);
I = zeros(length(t)+1, Nreps);
for k=1:length(t),
  C(k+1, :) = C(k, :) + (f_ext(k, :) - f_ext(k+1, :)).*exp(delta*(t_ext(k, :)+0.5));

  I(k+1, :) = I(k, :) + ...
      f_ext(k, :) / D .* (exp(D*(t_ext(k, :) + 0.5)) - exp(D*t_ext(k, :))) + ...
      C(k, :) / (D-delta) .* (exp((D-delta)*(t_ext(k, :) + 0.5)) - exp((D-delta)*t_ext(k, :))) + ...
      f_ext(k+1, :) / D .* (exp(D*t_ext(k+1, :)) - exp(D*(t_ext(k+1, :)-0.5))) + ...
      C(k+1, :) / (D-delta) .* (exp((D-delta)*t_ext(k+1, :)) - exp((D-delta)*(t_ext(k+1, :)-0.5)));
end
%[C, I]
%keyboard
m = mu + 1/delta * exp(-D* t_ext(2:end, :)) .* I(2:end, :);


function ll = gaussianLogLikelihood(x, y, vars),

ll = sum(-.5*numel(x)*log(2*pi) - .5*sum(log(vars)) ...
	 -.5*sum((x-y).^2 ./ vars));


function ll = ODEObjective(f, fvar, y, yvar, t, fest, params),

Nreps = size(f, 2);
if size(fest, 2) > 1,
  fest = fest';
end
m = integrateODE(fest, t, params);
ll = -gaussianLogLikelihood([m(:, ones(Nreps, 1)), fest(:, ones(Nreps, 1))], ...
			    [y, f], [yvar, fvar]);


function ll = ODEObjective2(f, fvar, y, yvar, t, fest, params),

m = integrateODE(fest, t, params);
ll = -gaussianLogLikelihood([m, fest], [y, f], [yvar, fvar]);
