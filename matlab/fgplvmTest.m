% FGPLVMTEST Test the gradients of the gpCovGrads function and the fgplvm models.

% FGPLVM

q = 2;
d = 3;
N = 10;
k = 5;
kernType = {'rbf', 'lin', 'rbfard', 'mlp', 'mlpard', 'white'};
kernType = {'rbf'};
Y = randn(N, d);
backType = 'mlp';
dynType = 'gp';
learn = 0; % dont' test learning of dynamics.
learnScales = 0; % don't test learning of dynamics scales.
diff = 1; % Use diffs for generating dynamics.

approxType = {'ftc', 'dtc', 'fitc', 'pitc'};
for back = [false true]
  for dyn = [false true];
    for a = 1:length(approxType)
      options = fgplvmOptions(approxType{a});
      options.kern = kernType;
      options.numActive = k;
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
        options.back = backType;
        options.backOptions = feval([backType 'Options']);
        options.optimiseInitBack = 0;
      end
      model = fgplvmCreate(q, d, Y, options);
      if dyn
        switch dynType 
         case 'gp'
          options = gpOptions(approxType{a});
          options.kern = kernType;
          options.numActive = k;
          model = fgplvmAddDynamics(model, 'gp', options, diff, learn);
         otherwise
           model = fgplvmAddDynamics(model, dynType);
           
        end
          
      end

      
      initParams = fgplvmExtractParam(model);
      % this creates some nasty parameters.
      initParams = randn(size(initParams))./randn(size(initParams));
      % This forces kernel computation.
      model = fgplvmExpandParam(model, initParams);
%      model.kern = kernSetWhite(model.kern, 1e-6);
%      params = fgplvmExtractParam(model);
%      model = fgplvmExpandParam(model, params);
      if dyn
        if strcmp(dynType, 'robOne')
          aveR = mean(model.dynamics.r);
          model.dynamics.b = model.dynamics.a/aveR;
        end
        if strcmp(dynType, 'robTwo')
          aveR = mean(model.dynamics.r);
          model.dynamics.b = model.dynamics.a/aveR;
        end
      end
      if ~back & ~dyn
        gpCovGradsTest(model);
      end
      gradientCheck(initParams, 'fgplvmObjective', 'fgplvmGradient', ...
               model);
    end
  end
end