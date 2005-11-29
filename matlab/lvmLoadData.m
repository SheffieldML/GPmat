function [Y, lbls] = lvmLoadData(dataset)

% LVMLOADDATA Load a dataset.

% DATASETS

lbls = [];
switch dataset
 
  %/~
 case 'run1'
  [Y, connect] = mocapLoadTextData('run1');
  Y = Y(1:4:end, :);
 case 'horse'
  load horse.dat;
  horse([133, 309], :) = [];
  ordinalAtr = [1 2 7 8 9 10 11 12 13 14 15 17 18 21];
  horse(:, ordinalAtr) = horse(:, ordinalAtr) - 1;

  % Set labels as outcome --- 1 lived, 2 died, 3 was euthanized
  lbls = zeros(size(horse, 1), 3);
  for i = 1:size(horse, 1)
    lbls(i, horse(i, 23)) = 1;
  end
  
  horse(:, 23:28) = [];
  horse(:, 3) = [];
  Y = horse;
  % Normalise Gaussian variables.
  for i = [3 4 5 15 18 19 21];
    va = var(Y(find(~isnan(Y(:, i))), i));
    Y(:, i) = Y(:, i)/va;
  end
  
  case 'pitch'
   load -ascii bballpitch.mat
   Y = bballpitch(:,8:end);

 case 'ratemaps'
  % Fix seeds
  Y = loadRateMap('ratemaps')';
  Y = Y(1:10000, :);
 
 case 'cepstral'
  Y = load('cepvecs');
  Y = Y(1:5000, :);
  %~/
 
 case 'brendan'
  load frey_rawface.mat
  Y = double(ff)';
  
 case 'digits'
  
  % Fix seeds
  randn('seed', 1e5);
  rand('seed', 1e5);

  load usps_train.mat
  % Extract 600 of digits 0 to 4
  [ALL_T, sortIndices] = sort(ALL_T);
  ALL_DATA = ALL_DATA(sortIndices(:), :);
  Y = [];
  lbls = [];
  numEachDigit = 600;
  for digit = 0:4;
    firstDigit = min(find(ALL_T==digit));
    Y = [Y; ALL_DATA(firstDigit:firstDigit+numEachDigit-1, :)];
    lbl = zeros(1, 5);
    lbl(digit+1) = 1;
    lbls = [lbls; repmat(lbl, numEachDigit, 1)];
  end

 case 'twos'  
  % load data
  load twos
  Y = 2*a-1;

 case 'oil'
  load 3Class.mat
  Y = DataTrn;
  lbls = DataTrnLbls;

 case 'oil100'
  randn('seed', 1e5);
  rand('seed', 1e5);
  load 3Class.mat
  Y = DataTrn;
  lbls = DataTrnLbls;
  indices = randperm(size(Y, 1));
  indices = indices(1:100);
  Y = Y(indices, :);
  lbls = lbls(indices, :);

 case 'swissRoll'
  load swiss_roll_data
  Y = X_data(:, 1:1000)';
 
end
