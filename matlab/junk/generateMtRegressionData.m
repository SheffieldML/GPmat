% GENERATEMTREGRESSIONDATA Tries to load a point-set regression sampled data set otherwise generates it.

% MTIVM

try
  load mtRegressionData.mat 
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
    if findstr(errMsg, 'read file')
      noFile = 1;
    end
  end
  if noFile
    fprintf('File not available ... generating data ... this will take some time\n');
    
    kernelType = 'rbf';
    numIn = 4;
    numTasks = 4;
    N = 2000;
    numDataPerTask = N*ones(1, numTasks);
    numTasksTest = 10;
    NTest = 500;
   numDataPerTaskTest = NTest*ones(1, numTasksTest);
   % Theta is of the form [rbfInverseWidth, rbfMultiplier, noiseVariance, biasVariance]
   trueTheta = [1 1 100 0];
   trueTheta = thetaConstrain(trueTheta);
   [X, y] = mtSampleData(trueTheta, kernelType, numIn, numTasks, numDataPerTask);
   [testX, testY] = mtSampleData(trueTheta, kernelType, numIn, numTasksTest, numDataPerTaskTest);

   save('mtRegressionData.mat', 'trueTheta', 'X', ...
	'y', 'testX', 'testY')
  else
    error(lasterr)
  end
end