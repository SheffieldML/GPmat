% VOWELMUTLIDEMO Try the multi-class IVM for classification of vowels.

prior = 0;
display = 1;
% Sample a regression data-set.
generateVowelData;
numTasks = size(vowelY, 2);
kernelType = 'rbf';
noiseType = 'multiprobit';
selectionCriterion = 'entropy';

for times = 1:10
  counter = 0;
  for d = [100 200 300 400 500 600 800 1000 1500 2000] %f2000
    counter = counter + 1;
    models(counter, times) = ivm(X, vowelY, ...
				 kernelType, noiseType, ...
				 selectionCriterion, d);
    
    
    % selection criteria
    for i = 1:4
      models(counter, times) = ivmOptimiseIVM(models(counter, times), display);
      models(counter, times) = ivmOptimiseKernel(models(counter, times), prior);
      disp(exp(models(counter, times).lntheta))
    end
    models(counter, times) = ivmOptimiseIVM(models(counter, times), display);
    out = zeros(size(vowelYTest));
    for taskNo = 1:numTasks
      [void, out(:, taskNo)] = ivmfwd(XTest, models(counter, times));
    end
    [void, classPred] = max(out, [], 2);
    [void, classTrue] = max(vowelYTest, [], 2);
    
    confusMat = zeros(11);
    for i = 1:length(classPred)
      confusMat(classPred(i), classTrue(i)) = confusMat(classPred(i), classTrue(i)) +1; 
    end
    
    perCentCorrect(counter, times) = sum(diag(confusMat))/sum(sum(confusMat))*100;
  end
end
