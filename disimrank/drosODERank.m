function [r, lls, fitparams] = drosODERank(drosexp, drosTF, tf, probes, t, doPlot),

% DROSODERANK Rank targets by fitting an ODE to raw expression data.
% FORMAT
% DESC Rank targets by fitting an ODE to raw expression data.
% ARG drosexp : the drosexp data structure from drosLoadData
% ARG drosTF : the drosTF data structure from drosLoadData
% ARG tf : TF symbol, should be in {'bap', 'bin', 'mef2', 'tin', 'twi'}
% ARG probes : probe sets to rank
% ARG t : time points to use (optional, default=all)
% ARG doPlot : Plot all models (optional, default=false)
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
  doPlot = 0;
end
N = length(drosexp.probes);

[tfprof, yvar, gene, times, scale, rawExp, rawVar] = drosGetData(drosexp, {drosTF.probes.(tf)});

tfprof = cat(2, tfprof{:});
tfprof = tfprof(t, :);

options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);

fitparams = {};
lls = zeros(size(probes));

for k=1:length(probes),
  fprintf('Scoring %s (%d/%d)...\n', probes{k}, k, length(probes));
  [y, yvar, gene, times, scale, rawExp, rawVar] = drosGetData(drosexp, probes(k));
  y = cat(2, y{:});
  yvar = cat(2, yvar{:});
  y = y(t, :);
  yvar = yvar(t, :);

  [fitparams{k}, lls(k)] = fminsearch(@(params) ODEObjective(tfprof, y, yvar, t, params), [0, 0, 0, 0], options);
  if doPlot,
    plot(t, tfprof, 'b'); hold on
    plot(t, y, 'r');
    plot(t, integrateODE(tfprof, t, fitparams{k}), 'g'); hold off;
    pause
  end
end

lls = -lls;
I = drosFindGeneinds(drosexp, probes, 0, 1);
[foo, J] = sort(lls, 'descend');
r = drosRemoveDuplicateGenes(drosexp, I(J));


function m = integrateODE(f, t, params),

Nreps = size(f, 2);
t_ext = [0, t];
t_ext = t_ext(ones(1, Nreps), :)';

B = exp(params(1));
D = exp(params(2));
S = exp(params(3));
delta = exp(params(4));

p = exp(-delta * t_ext) .* cumtrapz([zeros(1, Nreps); f] .* exp(delta * t_ext));
m = B/D + S * exp(-D * t_ext) .* cumtrapz(p .* exp(D * t_ext));
m = m(2:end, :);


function ll = gaussianLogLikelihood(x, y, vars),

ll = sum(-.5*numel(x)*log(2*pi) - .5*sum(log(vars)) ...
	 -.5*sum((x-y).^2 ./ vars));


function ll = ODEObjective(f, y, yvar, t, params),

m = integrateODE(f, t, params);
ll = -gaussianLogLikelihood(m, y, yvar);
