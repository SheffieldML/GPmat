function modelRet = gpTest

% GPTEST Test the gradients of the gpCovGrads function and the gp models.
% FORMAT
% DESC runs some tests on the GP-LVM code in the GP toolbox to
% test that it is working.
% RETURN model : a cell array of models used for testing.
%
% SEEALSO : modelTest
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2009


% GP

q = 2;
d = 3;
N = 10;
Nseq =4;
k = 5;
kernType = {'rbf', 'lin', 'rbfard', 'mlp', 'mlpard', 'white'};
kernType = 'rbf';
meanFunctionType = 'mlp';
learnScales = true; % test learning of output scales.
X = randn(N, q);
Yorig = randn(N, d);
indMissing = find(rand(N, d)>0.7);
%indMissing = [9 19 29];
approxType = {'ftc', 'dtc', 'dtcvar', 'fitc', 'pitc'};
approxType = {'dtcvar'};
counter = 0;
for optimiseBeta = [false true]
  for meanFunction = [false true]
    %  for missing = [false true]
    for missing = [false]
      for fixInducing = [false true] 
        Y = Yorig;
        if missing
          Y(indMissing) = NaN;
        end
        if meanFunction & missing
          continue
        end
        for a = 1:length(approxType)
          options = gpOptions(approxType{a});
          options.learnScales = learnScales;
          options.kern = kernType;
          options.numActive = k;
          options.isSpherical = ~missing;
          options.isMissingData = missing;
          options.fixInducing = fixInducing;
          options.optimiseBeta = optimiseBeta;
          if optimiseBeta & strcmp(approxType{a}, 'ftc')
            options.beta = 1000;
          end
          if meanFunction 
            disp(['Mean Function installed, ' ...
                  'with ' approxType{a} ...
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
          if ~optimiseBeta
            disp(['Beta not optimised.'])
          end
        
          if meanFunction
            options.meanFunction = meanFunctionType;
            options.meanFunctionOptions = feval([meanFunctionType 'Options']);
          end
          model = gpCreate(q, d, X, Y, options);
          
          
          initParams = gpExtractParam(model);
          % this creates some nasty parameters.
          initParams = randn(size(initParams));%./(100*randn(size(initParams)));
          
          % This forces kernel computation.
          model = gpExpandParam(model, initParams);
          gpCovGradsTest(model);
          modelGradientCheck(model);
          counter = counter + 1;
          modelRet{counter} = model;
        end
      end
    end
  end
end
