function [Y, lbls] = lvmLoadData(dataset)

% LVMLOADDATA Load a dataset.

% DATASETS

% get directory
baseDir = datasetsDirectory;

lbls = [];
switch dataset
 
 case 'robotWireless'
  Y = parseWirelessData([baseDir 'uw-floor.txt']);
  Y = Y(1:215, :);
 
 case 'robotWirelessTest'
  Y = parseWirelessData([baseDir 'uw-floor.txt']);
  Y = Y(216:end, :);

 case 'robotTraces'
  Y = csvread([baseDir 'Trace-3rdFloor-01.uwar.slam'], 1, 0); 
  Y = Y(1:floor(end/2), 4:end);
  Y(find(Y==-100))=-92;
  Y = (Y + 85)/15;
 case 'robotTracesTest'
  Y = csvread([baseDir 'Trace-3rdFloor-01.uwar.slam'], 1, 0); 
  Y = Y(ceil(end/2):end, 4:end);
  Y(find(Y==-100))=-92;
  Y = (Y + 85)/15;
  %/~
 case 'robotAtrium'
  Y = parseWirelessData([baseDir 'uw-atrium.txt']);
  Y = Y(1:281, :);
 
 case 'robotAtriumTest'
  Y = parseWirelessData([baseDir 'uw-atrium.txt']);
  Y = Y(282:end, :);

 case 'accWalkSitJog'
  try
    load([baseDir 'walkSitJog.mat']); 
  catch    
    [void, errid] = lasterr;
    if strcmp(errid, 'MATLAB:load:couldNotReadFile')
      skel = acclaimReadSkel([baseDir 'mocap\cmu\86\86.asf']);
      [channels, skel] = acclaimLoadChannels([baseDir 'mocap\cmu\86\86_10.amc'], skel);
      % Remove root z & x pos channels.
      channels(:, [1 3]) = zeros(size(channels, 1), 2);
      Y = channels(1:4:end, :);
      save([baseDir 'walkSitJog.mat'], 'Y')
    else
      error(lasterr);
    end
  end
 case 'horse'
  load([baseDir 'horse.dat']);
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
   load('-ascii', [baseDir 'bballpitch.mat']);
   Y = bballpitch(:,8:end);

 case 'ratemaps'
  % Fix seeds
  Y = loadRateMap([baseDir 'ratemaps'])';
  Y = Y(1:10000, :);
 
 case 'cepstral'
  Y = load([baseDir 'cepvecs']);
  Y = Y(1:5000, :);

 case 'vowels90'
  load([baseDir 'jon_vowel_data']);
  Y = [a_raw(1:30:end, :); ae_raw(1:30:end, :); ao_raw(1:30:end, :); ...
       e_raw(1:30:end, :); i_raw(1:30:end, :); ibar_raw(1:30:end, :); ...
       o_raw(1:30:end, :); schwa_raw(1:30:end, :); u_raw(1:30:end, :)];
  lbls = [];
  for i = 1:9
    lbl = zeros(1, 9);
    lbl(i) = 1;
    lbls = [lbls; repmat(lbl, 10, 1)];
  end
 case 'vowels3'
  load([baseDir 'jon_vowel_data']);
  Y = [a_raw(1:3:end, :); ae_raw(1:3:end, :); ao_raw(1:3:end, :); ...
       e_raw(1:3:end, :); i_raw(1:3:end, :); ibar_raw(1:3:end, :); ...
       o_raw(1:3:end, :); schwa_raw(1:3:end, :); u_raw(1:3:end, : ...
                                                    )];
  lbls = [];
  for i = 1:9
    lbl = zeros(1, 9);
    lbl(i) = 1;
    lbls = [lbls; repmat(lbl, 100, 1)];
  end
 case 'breakdance'
  Y = mocapLoadTextData([baseDir 'EricCamper04'], 0);
  Y = Y(1:4:end, :);
 
 case 'dance'
  Y = mocapLoadTextData([baseDir 'dancemotion'], 0);
  Y = Y(1:4:end, :);

 case 'danceReduced'
  Y = mocapLoadTextData([baseDir 'dancemotion'], 0);
  Y = Y(1:16:end, :);
  %~/
 case 'vowels'
  load([baseDir 'jon_vowel_data']);
  Y = [a_raw; ae_raw; ao_raw; ...
       e_raw; i_raw; ibar_raw; ...
       o_raw; schwa_raw; u_raw];
  Y(:, [13 26]) = [];
  lbls = [];
  for i = 1:9
    lbl = zeros(1, 9);
    lbl(i) = 1;
    lbls = [lbls; repmat(lbl, size(a_raw, 1), 1)];
  end
 

 case 'stick'
  Y = mocapLoadTextData([baseDir 'run1']);
  Y = Y(1:4:end, :);
 
 case 'walkSitJog'
  Y = mocapLoadTextData([baseDir 'walkSitJog']);
  Y = Y(1:4:end, :)*200;
 
 case 'brendan'
  load([baseDir 'frey_rawface.mat']);
  Y = double(ff)';
  
 case 'digits'
  
  % Fix seeds
  randn('seed', 1e5);
  rand('seed', 1e5);

  load([baseDir 'usps_train.mat']);
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
  load([baseDir 'twos']);
  Y = 2*a-1;

 case 'oil'
  load([baseDir '3Class.mat']);
  Y = DataTrn;
  lbls = DataTrnLbls;

 case 'oilTest'
  load([baseDir '3Class.mat']);
  Y = DataTst;
  lbls = DataTstLbls;

 case 'oilValid'
  load([baseDir '3Class.mat']);
  Y = DataVdn;
  lbls = DataVdnLbls;

 case 'oil100'
  randn('seed', 1e5);
  rand('seed', 1e5);
  load([baseDir '3Class.mat']);
  Y = DataTrn;
  lbls = DataTrnLbls;
  indices = randperm(size(Y, 1));
  indices = indices(1:100);
  Y = Y(indices, :);
  lbls = lbls(indices, :);

 case 'swissRoll'
  load([baseDir 'swiss_roll_data']);
  Y = X_data(:, 1:1000)';
 
 
 otherwise
  error('Unknown data set requested.')
 
end
