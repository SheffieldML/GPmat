function f = gpPosteriorSample(kernType, numSamps, params, lims, seed, bw)
  
% GPPOSTERIORSAMPLE Create a plot of samples from a posterior covariance.
% FORMAT
% DESC creates a plot of samples from a kernel with the given
% parameters and variance.
% ARG kernType : the type of kernel to sample from.
% ARG numSamps : the number of samples to take.
% ARG params : parameter vector for the kernel.
% ARG lims : limits of the x axis.
%
% COPYRIGHT : Neil D. Lawrence, 2008

% GP

global printDiagram
  
  
if nargin < 6
  bw = false;
  if nargin < 5
    seed = [];
    if nargin < 4
      lims = [-3 3];
      if nargin < 3
        params = [];
        if nargin < 2
          numSamps = 10;
        end
      end
    end
  end
end
if ~isempty(seed)
  randn('seed', seed);
  rand('seed', seed);
end
t_star = linspace(lims(1), lims(2), 200)';

kern = kernCreate(t_star, kernType);
if isfield(kern, 'comp')
  for i=1:length(kern.comp)
    kern.comp{i}.transforms = [];
  end
end

if ~isempty(params)
  feval = str2func([kern.type 'KernExpandParam']);
  kern = feval(kern, params);
end
feval = str2func([kern.type 'KernExtractParam']);
[params, names] = feval(kern);
paramStr = [];
for i = 1:length(names)
  Name = names{i};
  Name(1) = upper(Name(1));
  ind = find(Name==' ');
  Name(ind+1) = upper(Name(ind+1));
  Name(ind) = '';
  paramStr = [paramStr Name num2str(params(i))];
  
end
paramStr(find(paramStr==46)) = 'p';
infoStr = ['Samples' num2str(numSamps) 'Seed' num2str(randn('seed'))];

% Covariance of the prior.
K_starStar = kernCompute(kern, t_star, t_star);

% Generate "training data" from a sine wave.
t = rand(5, 1)*(lims(2)-lims(1))*0.6+lims(1)*0.6;
f = sin(t);

% Compute kernel for training data.
K_starf = kernCompute(kern, t_star, t);
K_ff = kernCompute(kern, t);

% Mean and covariance of posterior.
fbar = K_starf*pdinv(K_ff)*f;
Sigma = K_starStar - K_starf*pdinv(K_ff)*K_starf';

% Sample from the posterior.
fsamp = real(gsamp(fbar, Sigma, numSamps));

% Plot and save
figure
linHand = plot(t_star, fsamp);
hold on
linHandPlot = plot(t, f, 'r.')
set(linHandPlot, 'markersize', 30)

zeroAxes(gca, 0.025, 18, 'times');
set(linHand, 'linewidth', 1)
app = '';
if bw
  set(linHand, 'color', [0 0 0])
  set(linHandPlot, 'color', [0 0 0])
  app = 'bw';
end

if iscell(kernType)
  KernType = [];
  for i = length(kernType):-1:1
    KernType = [kernType{i} KernType];
    KernType(1) = upper(KernType(1));
  end
else
  KernType(1) = upper(kernType(1));
end

if exist('printDiagram', 'var') & printDiagram
  printPlot(['gpPosteriorSample' KernType infoStr paramStr app], ...
            '../tex/diagrams', '../html')
end


