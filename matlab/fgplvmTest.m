% FGPLVMTEST Test the gradients of the gpCovGrads function and the fgplvm models.

% FGPLVM

q = 2;
d = 3;
N = 10;
k = 5;
kernType = {'rbf', 'lin', 'rbfard', 'mlp', 'mlpard', 'white'};

Y = randn(N, d);


approxType = {'ftc', 'dtc', 'fitc', 'pitc'}; %{'dtc', 'fitc', 'pitc'};
for back = [false true]
  for dyn = [false true];
    for a = 1:length(approxType)
      if back & dyn
        disp(['Back constrained, ' ...
              'with dynamics and ' approxType{a} ...
              ' approximation.'])
      elseif dyn
        disp(['With dynamics and ' approxType{a} ...
              ' approximation.'])
      elseif back
        disp(['Back constrained and ' approxType{a} ...
              ' approximation.'])
      else
        disp([approxType{a} ' approximation.'])
        
      end
      if back
        model = fgplvmCreate(Y, q, approxType{a}, k, kernType);
      else
        model = fgplvmCreate(Y, q, approxType{a}, k, kernType, ...
                             'gaussian', 'linear');
      end
      if dyn
        model = fgplvmAddDynamics(model, kernType, 100);
      end
      
      initParams = fgplvmExtractParam(model);
      % this creates some nasty parameters.
      initParams = randn(size(initParams))./randn(size(initParams));
      % This forces kernel computation.
      model = fgplvmExpandParam(model, initParams);
      if ~back & ~dyn
        gpCovGradsTest(model);
      end
      gradientCheck(initParams, 'fgplvmObjective', 'fgplvmGradient', ...
               model);
    end
  end
end