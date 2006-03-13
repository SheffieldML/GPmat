% FGPLVMTEST Test the gradients of the gpCovGrads function and the fgplvm models.

% FGPLVM

% TODO: add checks for models with missing data.

q = 2;
d = 3;
N = 10;
k = 5;
kernType = {'rbf', 'lin', 'rbfard', 'mlp', 'mlpard', 'white'};
kernType = 'rbf'
backType = 'mlp';
dynType = 'gp';
learn = 0; % dont' test learning of dynamics.
learnScales = 0; % don't test learning of dynamics scales.
diff = 1; % Use diffs for generating dynamics.

Yorig = randn(N, d);
indMissing = find(rand(N, d)>0.9);
%Y(ind) = NaN;

approxType = {'ftc', 'dtc', 'fitc', 'pitc'};
approxType ={'fitc'}
for back = [false]% true]
  for missing = [false]% true]
    for fixInducing = [ false]% true] 
      Y = Yorig;
      if missing
        Y(indMissing) = NaN;
      end
      if back & missing
        continue
      end
      for dyn = [false]% true];
        for a = 1:length(approxType)
          options = fgplvmOptions(approxType{a});
          options.kern = kernType;
          options.numActive = k;
          optionsDyn = options;
          options.isSpherical = ~missing;
          options.isMissingData = missing;
          options.fixInducing = fixInducing;
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
          if missing
            disp(['Missing data used.'])
          end
          if fixInducing
            disp(['Inducing variables fixed.'])
            options.fixIndices = round(linspace(1, size(Y, 1), k));
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
              model = fgplvmAddDynamics(model, 'gp', optionsDyn, diff, learn);
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
          if ~back & ~dyn & ~missing
            gpCovGradsTest(model);
          end
          gradientCheck(initParams, 'fgplvmObjective', 'fgplvmGradient', ...
                        model);
        end
      end
    end
  end
end