% VIVMUSPSRESULTS Summarise the USPS result files in LaTeX.

% IVM

invariance = 'translate';
dVal = 1000;
try 
  load vivmUspsResults
catch
  [void, id] = lasterr;
  if strcmp(id, 'MATLAB:load:couldNotReadFile')
    [origX, origy, XTest, yTest] = mapLoadData('usps');
    [void, yTest] = max(yTest, [], 2);
    yTest = yTest - 1;
    
    for expNo = 1
      for seed = 1:10
        mu = zeros(size(yTest));
        for digit = 0:9
          load(['demUsps' num2str(digit) '_' invariance '_' num2str(expNo) '_d' ...
                num2str(dVal) '_seed' num2str(seed*1e5)]);
          [X, y] = ivmVirtual(origX(origIndex, :), origy(origIndex, :), invariance);
          
          testEr(seed, digit+1) = testError;
          model = ivmReconstruct(kern, noise, ivmInfo, ...
                                 X, y(:, digit+1));
          mu(:, digit+1) = ivmPosteriorMeanVar(model, XTest);
          mu(:, digit+1) = mu(:, digit+1) + model.noise.bias;
        end
        [void, yPred] = max(mu, [], 2);
        yPred = yPred - 1;
        overallError(seed) = 1 - sum(yPred == yTest)/size(yTest, 1);
        save('vivmUspsResults', 'testEr', 'overallError');
        
      end
    end
  else
    error(lasterr)
  end
end

% Now make the table
meanEr = squeeze(mean(testEr*100, 1));
varEr = squeeze(var(testEr*100, [], 1));
meanOverEr = mean(overallError*100);
varOverEr = var(overallError*100);
for i = 1:1
  for j = 1:10
    fprintf(['$%s \\pm %1.2f$'], numsf2str(meanEr(i, j), 3), sqrt(varEr(i, j)));
    fprintf(' & ');
  end
  fprintf(['$%s \\pm %1.2f$\\\n'], numsf2str(meanOverEr(i), 3), sqrt(varOverEr(i)))
end
