% IVMUSPSRESULTS Summarise the USPS result files in LaTeX.

% IVM

try 
  load ivmUspsResults
catch
  [void, id] = lasterr;
  if strcmp(id, 'MATLAB:load:couldNotReadFile')
    [X, y, XTest, yTest] = mapLoadData('usps');
    [void, yTest] = max(yTest, [], 2);
    yTest = yTest - 1;
    
    for expNo = 1:2
      for seed = 1:10
        mu = zeros(size(yTest));
        for digit = 0:9
          load(['demUsps' num2str(digit) '_' num2str(expNo) '_' ...
                num2str(seed*1e5)]);
          
          testEr(expNo, seed, digit+1) = testError;
          model = ivmReconstruct(kernStore, noiseStore, ivmInfoStore, ...
                                 X, y(:, digit+1));
          mu(:, digit+1) = ivmPosteriorMeanVar(model, XTest);
          mu(:, digit+1) = mu(:, digit+1) + model.noise.bias;
        end
        [void, yPred] = max(mu, [], 2);
        yPred = yPred - 1;
        overallError(expNo, seed) = 1 - sum(yPred == yTest)/size(yTest, 1);
        
      end
    end
    save('ivmUspsResults', 'testEr', 'overallError');
  else
    error(lasterr)
  end
end

% Now make the table
meanEr = squeeze(mean(testEr*100, 2));
varEr = squeeze(var(testEr*100, [], 2));
meanOverEr = mean(overallError*100, 2);
varOverEr = var(overallError*100, [], 2);
for i = 1:2
  for j = 1:10
    fprintf(['$%s \\pm %1.2f$'], numsf2str(meanEr(i, j), 3), sqrt(varEr(i, j)));
    fprintf(' & ');
  end
  fprintf(['$%s \\pm %1.2f$\n'], numsf2str(meanOverEr(i), 3), sqrt(varOverEr(i)))
end
