% GTMOIL For visualising oil data --- uses NETLAB toolbox.


rand('state', 1e5);
randn('state', 1e5);

load 3Class

dataDim = size(DataTrn, 2);
latentDim = 2;


latentGridDims = [15 15]; 
numLatent = prod(latentGridDims);  % Number of latent points
numCentres = 16;


% Create and initialise GTM model
model = gtm(latentDim, numLatent, dataDim, numCentres, ...
   'gaussian', 0.1);

options = foptions;
options(7) = 1;   
model = gtminit(model, options, DataTrn, 'regular', latentGridDims, [4 4]);

options = foptions;
options(14) = 1000;
options(1) = 1;

[model, options] = gtmem(model, DataTrn, options);

% Plot posterior means
X = gtmlmean(model, DataTrn);
figure, hold on

for i = 1:size(X, 1)
  if(DataTrnLbls(i, 1) == 1)
    plot(X(i, 1), X(i, 2), 'r+')
  elseif(DataTrnLbls(i, 2) == 1)
    plot(X(i, 1), X(i, 2), 'bo')
  else
    plot(X(i, 1), X(i, 2), 'mx')
  end
end
