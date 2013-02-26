function [X, y, XTest, yTest] = ncnmLoadData(dataSet,seedVal)

% NCNMLOADDATA Load a dataset.
%
%	Description:
%
%	[X, Y] = NCNMLOADDATA(DATASETNAME) loads particular data set for use
%	with NCNM demos.
%	 Returns:
%	  X - input values for data.
%	  Y - target values for data.
%	 Arguments:
%	  DATASETNAME - the name of the data set to load.
%
%	[X, Y] = NCNMLOADDATA(DATASETNAME, SEEDVAL) loads particular data
%	set for use with NCNM demos initialising with a given random seed.
%	 Returns:
%	  X - input values for data.
%	  Y - target values for data.
%	 Arguments:
%	  DATASETNAME - the name of the data set to load.
%	  SEEDVAL - seed value to initialise with when generating data.
%
%	[X, Y, XTEST, YTEST] = NCNMLOADDATA(DATASETNAME) loads particular
%	data set for use with NCNM demos, this version will return test data
%	as well as training data.
%	 Returns:
%	  X - input values for data.
%	  Y - target values for data.
%	  XTEST - test input values for data.
%	  YTEST - test target values for data.
%	 Arguments:
%	  DATASETNAME - the name of the data set to load.
%	
%
%	See also
%	MAPLOADDATA, DEMTHREEFIVE, DEMUNLABELLED1, DEMPROBIT1


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


if nargin < 2
  seedVal = 1e5;
end
randn('seed', seedVal)
rand('seed', seedVal)
XTest = [];
yTest = [];
[dataSetStub, prob] = strtok(dataSet, '_');
prob = str2num(prob(2:end));
switch dataSetStub
  case 'threeFive'
   
   load usps_train
   X = ALL_DATA;
   y = ALL_T;
   load usps_test
   XTest = ALL_DATA;
   yTest = ALL_T;
   classTrue = 3;
   for i = [0 1 2 4 6 7 8 9];
     index = find(y == i);
     X(index, :) = [];
     y(index, :) = [];
     index = find(yTest == i);
     XTest(index, :) = [];
     yTest(index, :) = [];
   end
   y = (y == classTrue)*2 - 1;
   yTest = (yTest == classTrue)*2 - 1;
  case 'fourNine'
   
   load usps_train
   X = ALL_DATA;
   y = ALL_T;
   load usps_test
   XTest = ALL_DATA;
   yTest = ALL_T;
   classTrue = 4;
   for i = [0 1 2 3 5 6 7 8];
     index = find(y == i);
     X(index, :) = [];
     y(index, :) = [];
     index = find(yTest == i);
     XTest(index, :) = [];
     yTest(index, :) = [];
   end
   y = (y == classTrue)*2 - 1;
   yTest = (yTest == classTrue)*2 - 1;
 case 'thorsten'
  [y, X] = svmlread('example2/train_transduction.dat');
  [yTest, XTest] = svmlread('example2/test.dat');

 case 'usps'
  load usps_train
  X = ALL_DATA;
  range =  min(ALL_T):max(ALL_T);
  for i = 1:length(range)
    y(:, i) = (ALL_T == range(i))*2 - 1;
  end
  if nargout > 2
    load usps_test
    XTest = ALL_DATA;
    range =  min(ALL_T):max(ALL_T);
    for i = 1:length(range)
      yTest(:, i) = (ALL_T == range(i))*2 - 1;
    end
  end
  indUnlabelled = find(rand(size(y, 1)>prob));
  y(indUnlabelled, :) = NaN;
 case {'usps0', 'usps1', 'usps2', 'usps3', 'usps4', 'usps5', 'usps6', 'usps7', 'usps8', 'usps9'}
  digitNo = str2num(dataSetStub(end));
  load usps_train
  X = ALL_DATA;
  range =  min(ALL_T):max(ALL_T);
  for i = 1:length(range)
    y(:, i) = (ALL_T == range(i))*2 - 1;
  end
  if nargout > 2
    load usps_test
    XTest = ALL_DATA;
    range =  min(ALL_T):max(ALL_T);
    for i = 1:length(range)
      yTest(:, i) = (ALL_T == range(i))*2 - 1;
    end
  end
  y = y(:, digitNo+1);
  indUnlabelled = find(rand(size(y, 1), 1)>prob);
  y(indUnlabelled, :) = NaN;
  if nargout>2
    yTest = yTest(:, digitNo+1);
  end
       
end
