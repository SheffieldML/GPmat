% GENERATECLASSIFICATIONDATA Tries to load a sampled data set otherwise generates it.
 
seed = 1e5;
try
  load classificationData.mat 
catch

  noFile = 0;
  verString = version;

  if str2double(verString(1:3)) > 6.1
    [void, errid] = lasterr;
    if strcmp(errid, '')
      noFile = 1;
    end
  else
    errMsg = lasterr;
    if findstr(errMsg, 'file does not exist')
      noFile = 1;
    end
  end
  if noFile
    randn('seed', seed)
    rand('seed', seed)
    numIn = 2;
    N = 500;
    
    X = zeros(N, numIn);
    X = rand(N, numIn);

    kern = kernelCreate(X, 'rbf');
    kern.variance = 100;
    kern.inverseWidth = 10;
    
    K = kernCompute(kern, X);
    u = real(gaussSamp(K, 1)');
  
    p = cumGaussian(u);
    y = 2*(rand(size(u))>p)-1;
  end

  save('classificationData.mat', 'numIn', 'N', 'X', 'u', 'y')
end
