function [X, y, XTest, yTest] = ivmLoadData(dataset)

% IVMLOADDATA Load a dataset.

% IVM

switch dataset
 case 'usps'
  load ../data/usps_train
  X = ALL_DATA;
  range =  min(ALL_T):max(ALL_T);
  for i = 1:length(range)
    y(:, i) = (ALL_T == range(i))*2 - 1;
  end
  if nargout > 2
    load ../data/usps_test
    XTest = ALL_DATA;
    range =  min(ALL_T):max(ALL_T);
    for i = 1:length(range)
      yTest(:, i) = (ALL_T == range(i))*2 - 1;
    end
  end

end
