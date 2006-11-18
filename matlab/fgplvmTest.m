function modelRet = fgplvmTest

% FGPLVMTEST Test the gradients of the gpCovGrads function and the fgplvm models.
% FORMAT
% DESC runs some tests on the GP-LVM code in the FGPLVM toolbox to
% test that it is working.
% RETURN model : a cell array of models used for testing.
%
% SEEALSO : modelTest
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006


% FGPLVM

q = 2;
d = 3;
N = 10;
Nseq =4;
k = 5;
kernType = {'rbf', 'lin', 'rbfard', 'mlp', 'mlpard', 'white'};
kernType = 'rbf';
backType = 'mlp';
dynType = 'gpTime';
learn = true; % dont' test learning of dynamics.
diff = false; % Use diffs for generating dynamics.
seq(1) = 5;
seq(2) = 10;
learnScales = true; % test learning of output scales.

Yorig = randn(N, d);
indMissing = find(rand(N, d)>0.7);
%indMissing = [9 19 29];
approxType = {'ftc', 'dtc', 'fitc', 'pitc'};
counter = 0;
for back = [false true]
%  for missing = [false true]
  for missing = [false]
    for fixInducing = [false true] 
      Y = Yorig;
      if missing
        Y(indMissing) = NaN;
      end
      if back & missing
        continue
      end
      for dyn = [false true];
        for a = 1:length(approxType)
          options = fgplvmOptions(approxType{a});
          options.learnScales = learnScales;
          options.kern = kernType;
          options.numActive = k;
          options.isSpherical = ~missing;
          options.isMissingData = missing;
          optionsDyn = options;
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
              model = fgplvmAddDynamics(model, 'gp', optionsDyn, ...
                                        diff, learn, seq);
             case 'gpTime'
              t = [1:size(Y, 1)]';
              model = fgplvmAddDynamics(model, 'gpTime', optionsDyn, ...
                                        t, diff, learn, seq);
             otherwise
              model = fgplvmAddDynamics(model, dynType);
              
            end
            
          end
          
          
          initParams = fgplvmExtractParam(model);
          % this creates some nasty parameters.
          initParams = randn(size(initParams));%./(100*randn(size(initParams)));
          % This forces kernel computation.
          model = fgplvmExpandParam(model, initParams);
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
          seqX = randn(Nseq, model.q);
          seqY = randn(Nseq, model.d);
          seqY(find(rand(Nseq, model.d)>0.7)) = NaN;
          gradientCheck(seqX(:)', 'fgplvmSequenceObjective', 'fgplvmSequenceGradient', ...
                        model, seqY);
          counter = counter + 1;
          modelRet{counter} = model;
        end
      end
    end
  end
end