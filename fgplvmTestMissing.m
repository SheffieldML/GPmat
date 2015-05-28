% FGPLVMTESTMISSING Make sure missing data likelihood match full ones.

% FGPLVM

q = 2;
d = 3;
N = 10;
k = 5;
kernType = {'rbf', 'lin', 'rbfard', 'mlp', 'mlpard', 'white'};

Y = randn(N, d);

approxType = {'ftc', 'dtc', 'fitc', 'pitc'};
counter = 0;
for a = 1:length(approxType)
  for missing = [false true]
    counter = counter + 1;
    options = fgplvmOptions(approxType{a});
    options.kern = kernType;
    options.numActive = k;
    options.isSpherical = ~missing;
    options.isMissingData = missing;
    disp([approxType{a} ' approximation.'])
    if missing
      disp(['Missing data used.'])
    end
    model{counter} = fgplvmCreate(q, d, Y, options);
          
    if ~missing
      % Set up parameters for each model{counter}
      initParams = fgplvmExtractParam(model{counter});
      % this creates some nasty parameters.
      initParams = randn(size(initParams))./randn(size(initParams));
    end
    model{counter} = fgplvmExpandParam(model{counter}, initParams);
    % This forces kernel computation.
    ll = gpLogLikelihood(model{counter});
    disp(['Log likelihood ' num2str(ll)])
  end
end
