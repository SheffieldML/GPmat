funcs={'ngaussian', 'cumGaussian', 'invCumGaussian', 'gradLnCumGaussian', ...
       'lnCumGaussian', 'sigmoid', 'invSigmoid', 'erfcinv', 'gamma', ...
       'gammaln', 'digamma'}
for i=1:length(funcs)
  x = randn(1); 
  switch(funcs{i})
   case 'invSigmoid'
    x = rand(1);
   case 'invCumGaussian'
    x = rand(1);
   case 'erfcinv'
    x = rand(1);
   case 'gamma'
    x = abs(x);
   case 'gammaln'
    x = abs(x);
  end
  switch(funcs{i})
   case 'gradLnCumGaussian'
    y = gradLogCumGaussian(x);
   otherwise 
    y = feval(funcs{i}, x);
  end
  save([funcs{i} 'Test.mat'], 'x', 'y');
end
