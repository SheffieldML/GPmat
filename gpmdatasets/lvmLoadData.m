function [Y, lbls, Ytest, lblstest] = lvmLoadData(dataset, seedVal)

% LVMLOADDATA Load a latent variable model dataset.
% FORMAT
% DESC loads a data set for a latent variable modelling problem.
% ARG dataset : the name of the data set to be loaded. Currently
% the possible names are 'robotWireless',
% 'robotTwoLoops', 'robotTraces', 'robotTracesTest', 'cmu35gplvm',
% 'cmu35Taylor', 'cmu35walkJog', 'vowels', 'stick', 'brendan',
% 'digits', 'twos', 'oil', 'oilTest', 'oilValid', 'oil100',
% 'swissRoll', 'missa'.
% RETURN Y : the training data loaded in.
% RETURN lbls : a set of labels for the data (if there are no
% labels it is empty).
% RETURN Ytest : the test data loaded in. If no test set is
% available it is empty.
% RETURN lblstest : a set of labels for the test data (if there are
% no labels it is empty).
%
% SEEALSO : mapLoadData, datasetsDirectory
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2008, 2009

% DATASETS

  if nargin > 1
    randn('seed', seedVal)
    rand('seed', seedVal)
  end

  % get directory

  baseDir = datasetsDirectory;
  dirSep = filesep;
  lbls = [];
  lblstest = [];
  switch dataset
   case 'YaleSubset6_1'
          load([baseDir 'YaleSubset6.mat']);
          lbls = [height;width];
   case 'missa'
       load([baseDir 'miss-americaHD.mat']);
       lbls = [288 ;360];
   case 'movielens'
    try 
      load([baseDir 'movielens.mat']);
    catch
      [void, errid] = lasterr;
      if strcmp(errid, 'MATLAB:load:couldNotReadFile');
        fileName = [baseDir dirSep 'movielens' dirSep 'large' dirSep 'ratings.dat'];
        [users, films, ratings, timeStamp] = textread(fileName, '%n::%n::%n::%n');
        ind = randperm(size(users, 1));
        users = users(ind, :);
        films = films(ind, :);
        ratings = ratings(ind, :);
        numUsers = max(users);
        numFilms = max(films);
        
        numRatings = size(users, 1);
        numTrainRatings = ceil(0.8*numRatings);
        Y = spalloc(numFilms, numUsers, numTrainRatings);
        Ytest = spalloc(numFilms, numUsers, numRatings-numTrainRatings);
        indTrain = sub2ind(size(Y), films(1:numTrainRatings), users(1:numTrainRatings));
        indTest = sub2ind(size(Ytest), films(numTrainRatings+1:numRatings), users(numTrainRatings+1:numRatings));
        Y(indTrain) = ratings(1:numTrainRatings);
        Ytest(indTest) = ratings(numTrainRatings+1:numRatings);
        save([baseDir 'movielens.mat'], 'Y', 'lbls', 'Ytest', 'lblstest');
      else
        error(lasterr);
      end
    end
      
   case {'movielensSmall1', 'movielensSmall2', 'movielensSmall3', 'movielensSmall4', 'movielensSmall5'}

    partNo = str2num(dataset(end));
    uTrain = load([baseDir dirSep 'movielens' dirSep 'small' dirSep 'u' num2str(partNo) '.base']);
    numUsers = max(uTrain(:, 1));
    numFilms = max(uTrain(:, 2));
    numRatings = size(uTrain, 1);
    Y = spalloc(numFilms, numUsers, numRatings);
    
    for i = 1:size(uTrain, 1);
      Y(uTrain(i, 2), uTrain(i, 1)) = uTrain(i, 3);
    end
    meanY = mean(Y(find(Y)));
    Y(find(Y)) = (Y(find(Y))-meanY);
    uTest = load([baseDir dirSep 'movielens' dirSep 'small' dirSep 'u' num2str(partNo) '.test']);
    numTestRatings = size(uTest, 1);
    Ytest = spalloc(numFilms, numUsers, numTestRatings);
    for i = 1:size(uTest, 1);
      Ytest(uTest(i, 2), uTest(i, 1)) = uTest(i, 3);
    end
    Ytest(find(Ytest)) = (Ytest(find(Ytest))-meanY);
    stdY = sqrt(var(Y(find(Y))));
    Y = Y/stdY;
    Ytest = Ytest/stdY;
    
   case 'crabs'
    fid = fopen([baseDir 'crabs.dat']);
    if fid == -1
      error(['No such file name ' fileName])
    end
    readLine = fgets(fid);
    counter = 0;
    data = [];
    while readLine ~= -1
      readLine = fgets(fid);
      counter = counter + 1;
      parts = stringSplit(readLine, ' ');
      if parts{1} == 'B' && parts{2} == 'F'
        lbls(counter, :) = [0 0 0 1];
      elseif parts{1} == 'B' && parts{2} == 'M'
        lbls(counter, :) = [0 0 1 0];
      elseif parts{1} == 'O' && parts{2} == 'F'
        lbls(counter, :) = [0 1 0 0];
      else
        lbls(counter, :) = [1 0 0 0];
      end
      counter2 = 0;
      for i = 3:size(parts, 2)
        if ~isempty(str2num(parts{i}))
          counter2 = counter2 + 1;
          data(counter, counter2) = str2num(parts{i});
        end
      end
    end
    fclose(fid);
    data = data(:, 2:end);
    Y = data./repmat(sum(data, 2), 1, size(data, 2));
    
   case 'spellman'
    try 
      load([baseDir 'spellman.mat']);
      
    catch
      [void, errid] = lasterr;
      if strcmp(errid, 'MATLAB:load:couldNotReadFile');
        fileName = ([baseDir 'combined.csv']);
        [geneName, ...
         y1, y2, y3, y4, y5, y6, y7, y8, y9, ...
         y10, y11, y12, y13, y14, y15, y16, y17, y18, ...
         y19, y20, y21, y22, y23, y24, y25, y26, y27, ...
         y28, y29, y30, y31, y32, y33, y34, y35, y36, ...
         y37, y38, y39, y40, y41, y42, y43, y44, y45, ...
         y46, y47, y48, y49, y50, y51, y52, y53, y54, ...
         y55, y56, y57, y58, y59, y60, y61, y62, y63, ...
         y64, y65, y66, y67, y68, y69, y70, y71, y72, ...
         y73, y74, y75, y76, y77] = ...
            textread(fileName, ['%q ' ...
                            '%f %f %f %f %f %f %f %f %f' ...
                            '%f %f %f %f %f %f %f %f %f' ...
                            '%f %f %f %f %f %f %f %f %f' ...
                            '%f %f %f %f %f %f %f %f %f' ...
                            '%f %f %f %f %f %f %f %f %f' ...
                            '%f %f %f %f %f %f %f %f %f' ...
                            '%f %f %f %f %f %f %f %f %f' ...
                            '%f %f %f %f %f %f %f %f %f' ...
                            '%f %f %f %f %f'], ...
                     'headerlines', 1, 'whitespace', '', 'delimiter', ',');
        
        Y = [y1 y2 y3 y4 y5 y6 y7 y8 y9 ...
             y10, y11, y12, y13, y14, y15, y16, y17, y18, ...
             y19, y20, y21, y22, y23, y24, y25, y26, y27, ...
             y28, y29, y30, y31, y32, y33, y34, y35, y36, ...
             y37, y38, y39, y40, y41, y42, y43, y44, y45, ...
             y46, y47, y48, y49, y50, y51, y52, y53, y54, ...
             y55, y56, y57, y58, y59, y60, y61, y62, y63, ...
             y64, y65, y66, y67, y68, y69, y70, y71, y72, ...
             y73, y74, y75, y76, y77]';
        lbls = [repmat([1 0 0 0 0], 4, 1); 
                repmat([0 1 0 0 0], 18,1);
                repmat([0 0 1 0 0], 24,1);
                repmat([0 0 0 1 0], 17,1);
                repmat([0 0 0 0 1], 14,1)];
        save([baseDir 'spellman.mat'], 'Y', 'lbls');
      else
        error(lasterr);
      end
    end
   case 'robotWireless'
    Ydat = parseWirelessData([baseDir 'uw-floor.txt']);
    Ydat = (Ydat + 85)/15;
    lbls = 'connect';
    Y = Ydat(1:215, :);    
    Ytest = Ydat(216:end, :);

   case 'robotTwoLoops'
    Y = csvread([baseDir 'TwoLoops.slam'], 1, 0);
    Y = Y(1:floor(end/2), 4:end);
    Y(find(Y==-100))=-NaN;
    Y = (Y + 85)/15;
    lbls = 'connect';
    
   case 'robotTraces'
    Y = csvread([baseDir 'Trace-3rdFloor-01.uwar.slam'], 1, 0); 
    Y = Y(1:floor(end/2), 4:end);
    Y(:, [3 4 38]) = []; % Remove columns of missing data.
    Y(find(Y==-100))=NaN;
    Y = (Y + 85)/15;
    lbls = 'connect';

   case 'robotTracesTest'
    Y = csvread([baseDir 'Trace-3rdFloor-01.uwar.slam'], 1, 0); 
    Y = Y(ceil(end/2):end, 4:end);
    Y(:, [3 4 38]) = []; % Remove columns of missing data.
    Y(find(Y==-100))=NaN;
    Y = (Y + 85)/15;
    lbls = 'connect';

   case 'osr'
    t = linspace(-pi, pi, 256);
    window = (cos(t) + 1)/2;
    [x, fs] = wavread([baseDir 'osr' dirSep 'OSR_uk_000_0020_8k.wav']);
    [r, theta] = stft(x, fs, 256, window);
    Y = [r(:, 1:32) theta(:, 1:32)];

   case 'cmu49BalanceArm'
    [Y, lbls, Ytest, lblstest] = lvmLoadData('cmu49Balance');
    % Extract left arm only 
    ind = [41:47 49:50];
    Y = Y(:, ind);  
    Ytest = Ytest(:, ind);
    
   case 'cmu49Balance'
    
    try 
      load([baseDir 'cmu49Balance.mat']);
    catch
      [void, errid] = lasterr;
      if strcmp(errid, 'MATLAB:load:couldNotReadFile');
        sampleEvery = 32;
        skel = acclaimReadSkel([baseDir 'mocap' dirSep 'cmu' dirSep '49' dirSep '49.asf']);
        examples = {'18', '19'};
        testExamples = {'20'};
        % Label differently for each sequence
        exlbls = eye(length(examples));
        testexlbls = eye(length(testExamples));
        totLength = 0;
        totTestLength = 0;
        for i = 1:length(examples)
          [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep '49' dirSep '49_' ...
                              examples{i} '.amc'], skel);
          tY{i} = tmpchan(1:sampleEvery:end, :);
          tlbls{i} = repmat(exlbls(i, :), size(tY{i}, 1), 1);
          totLength = totLength + size(tY{i}, 1);
        end
        Y = zeros(totLength, size(tY{1}, 2));
        lbls = zeros(totLength, size(tlbls{1}, 2));
        endInd = 0;
        for i = 1:length(tY)
          startInd = endInd + 1;
          endInd = endInd + size(tY{i}, 1);
          Y(startInd:endInd, :) = tY{i};
          lbls(startInd:endInd, :) = tlbls{i};
        end
        for i = 1:length(testExamples)
          [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep '49' dirSep '49_' ...
                              testExamples{i} '.amc'], skel);
          tYtest{i} = tmpchan(1:sampleEvery:end, :);
          tlblstest{i} = repmat(testexlbls(i, :), size(tYtest{i}, 1), 1);
          totTestLength = totTestLength + size(tYtest{i}, 1);
        end
        Ytest = zeros(totTestLength, size(tYtest{1}, 2));
        lblstest = zeros(totTestLength, size(tlblstest{1}, 2));
        endInd = 0;
        for i = 1:length(tYtest)
          startInd = endInd + 1;
          endInd = endInd + size(tYtest{i}, 1);
          Ytest(startInd:endInd, :) = tYtest{i};
          lblstest(startInd:endInd, :) = tlblstest{i};
        end
        save([baseDir 'cmu49Balance.mat'], 'Y', 'lbls', 'Ytest', 'lblstest');
      else
        error(lasterr);
      end
    end
   
   case 'cmu07WalkFeet'
    [Y, lbls, Ytest, lblstest] = lvmLoadData('cmu07Walk');
    % Extract left arm only 
    ind = [7:12 14:19];
    Y = Y(:, ind);  
    Ytest = Ytest(:, ind);
    
   case 'cmu07Walk'
    
    try 
      load([baseDir 'cmu07Walk.mat']);
    catch
      [void, errid] = lasterr;
      if strcmp(errid, 'MATLAB:load:couldNotReadFile');
        sampleEvery = 4;
        skel = acclaimReadSkel([baseDir 'mocap' dirSep 'cmu' dirSep '07' dirSep '07.asf']);
        examples = {'01', '02'};
        testExamples = {'03'};
        % Label differently for each sequence
        exlbls = eye(length(examples));
        testexlbls = eye(length(testExamples));
        totLength = 0;
        totTestLength = 0;
        for i = 1:length(examples)
          [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep '07' dirSep '07_' ...
                              examples{i} '.amc'], skel);
          tY{i} = tmpchan(1:sampleEvery:end, :);
          tlbls{i} = repmat(exlbls(i, :), size(tY{i}, 1), 1);
          totLength = totLength + size(tY{i}, 1);
        end
        Y = zeros(totLength, size(tY{1}, 2));
        lbls = zeros(totLength, size(tlbls{1}, 2));
        endInd = 0;
        for i = 1:length(tY)
          startInd = endInd + 1;
          endInd = endInd + size(tY{i}, 1);
          Y(startInd:endInd, :) = tY{i};
          lbls(startInd:endInd, :) = tlbls{i};
        end
        for i = 1:length(testExamples)
          [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep '07' dirSep '07_' ...
                              testExamples{i} '.amc'], skel);
          tYtest{i} = tmpchan(1:sampleEvery:end, :);
          tlblstest{i} = repmat(testexlbls(i, :), size(tYtest{i}, 1), 1);
          totTestLength = totTestLength + size(tYtest{i}, 1);
        end
        Ytest = zeros(totTestLength, size(tYtest{1}, 2));
        lblstest = zeros(totTestLength, size(tlblstest{1}, 2));
        endInd = 0;
        for i = 1:length(tYtest)
          startInd = endInd + 1;
          endInd = endInd + size(tYtest{i}, 1);
          Ytest(startInd:endInd, :) = tYtest{i};
          lblstest(startInd:endInd, :) = tlblstest{i};
        end
        save([baseDir 'cmu07Walk.mat'], 'Y', 'lbls', 'Ytest', 'lblstest');
      else
        error(lasterr);
      end
    end 
    
    
   case 'cmu35gplvm'
    [Y, lbls, Ytest, lblstest] = lvmLoadData('cmu35WalkJog');
    skel = acclaimReadSkel([baseDir 'mocap' dirSep 'cmu' dirSep '35' dirSep '35.asf']);
    [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep '35' dirSep '35_01.amc'], skel);

    Ytest = Ytest(find(lblstest(:, 2)), :);
    lblstest = lblstest(find(lblstest(:, 2)), 2);

    %left indices
    xyzInd = [2];
    xyzDiffInd = [1 3];
    rotInd = [4 6];
    rotDiffInd = [5];
    generalInd = [7:38 41:47 49:50 53:59 61:62];

    jointAngles  = asin(sin(pi*Y(:, generalInd)/180));
    jointAnglesTest  = asin(sin(pi*Ytest(:, generalInd)/180));
    
    endInd = [];
    for i = 1:size(lbls, 2)
      endInd = [endInd max(find(lbls(:, i)))];
    end
    catJointAngles = [];
    xyzDiff = [];
    catSinCos = [];
    startInd = 1;
    for i = 1:length(endInd)
      ind1 = startInd:endInd(i)-1;
      ind2 = startInd+1:endInd(i);
      catJointAngles = [catJointAngles; ...
                        jointAngles(ind2, :)];
      xyzDiff = [xyzDiff;
                 Y(ind1, xyzDiffInd) - Y(ind2, xyzDiffInd) ...
                 Y(ind2, xyzInd)];
      catSinCos = [catSinCos; ...
                   sin(pi*Y(ind2, rotInd)/180) ...
                   sin(pi*Y(ind1, rotDiffInd)/180)-sin(pi*Y(ind2, rotDiffInd)/180) ...
                   cos(pi*Y(ind2, rotInd)/180) ...
                   cos(pi*Y(ind1, rotDiffInd)/180)-cos(pi*Y(ind2, rotDiffInd)/180)];
      startInd = endInd(i)+1;
    end
    Y = [catJointAngles xyzDiff catSinCos];
    lbls = [];
    
    endInd = [];
    for i = 1:size(lblstest, 2)
      endInd = [endInd max(find(lblstest(:, i)))];
    end
    catJointAnglesTest = [];
    xyzDiffTest = [];
    catSinCosTest = [];
    startInd = 1;
    for i = 1:length(endInd)
      ind1 = startInd:endInd(i)-1;
      ind2 = startInd+1:endInd(i);
      catJointAnglesTest = [catJointAnglesTest; ...
                          jointAnglesTest(ind2, :)];
      xyzDiffTest = [xyzDiffTest;
                     Ytest(ind1, xyzDiffInd) - Ytest(ind2, xyzDiffInd) ...
                     Ytest(ind2, xyzInd)];
      catSinCosTest = [catSinCosTest; ...
                       sin(pi*Ytest(ind2, rotInd)/180) ...
                       sin(pi*Ytest(ind1, rotDiffInd)/180)-sin(pi*Ytest(ind2, rotDiffInd)/180) ...
                       cos(pi*Ytest(ind2, rotInd)/180) ...
                       cos(pi*Ytest(ind1, rotDiffInd)/180)-cos(pi*Ytest(ind2, rotDiffInd)/180)];
      startInd = endInd(i)+1;
    end                                                
    Ytest = [catJointAnglesTest xyzDiffTest catSinCosTest];
    lblstest = [];

   case 'cmu35Taylor'
    % The CMU 35 data set as Graham Taylor et al. have it in their
    % NIPS 2006 paper.
    [rawY, lbls, rawYtest, lblstest] = lvmLoadData('cmu35WalkJog');
    rawYtest = rawYtest(45:end, :);
    lblstest = lblstest(45:end, 2);
    skel = acclaimReadSkel([baseDir 'mocap' dirSep 'cmu' dirSep '35' dirSep '35.asf']);
    [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep '35' dirSep '35_01.amc'], skel);
    posInd = [1 2 3];
    rotInd = [4 5 6];
    generalInd = [7:38 41:47 49:50 53:59 61:62];
    
    rawY = [rawY(:, posInd) sind(rawY(:, rotInd)) cosd(rawY(:, rotInd)) asind(sind(rawY(:, generalInd)))];
    rawYtest = [rawYtest(:, posInd) sind(rawYtest(:, rotInd)) cosd(rawYtest(:, rotInd)) asind(sind(rawYtest(:, generalInd)))];
    
    endInd = [];
    for i = 1:size(lbls, 2)
      endInd = [endInd max(find(lbls(:, i)))];
    end
    Y = [];
    startInd = 1;
    for i = 1:length(endInd)
      ind1 = startInd:endInd(i)-1;
      ind2 = startInd+1:endInd(i);
      Y = [Y; rawY(ind1, :) rawY(ind2, :)];
      startInd = endInd(i)+1;
    end  
    lbls = [];
    
    endInd = [];
    for i = 1:size(lblstest, 2)
      endInd = [endInd max(find(lblstest(:, i)))];
    end
    Ytest = [];
    startInd = 1;
    for i = 1:length(endInd)
      ind1 = startInd:endInd(i)-1;
      ind2 = startInd+1:endInd(i);
      Ytest = [Ytest; rawYtest(ind1, :) rawYtest(ind2, :)];
      startInd = endInd(i)+1;
    end                                                
    lblstest = [];

   case {'dur', 'cmp'}
    [Y, void, lbls] = synthLoadData(dataset);
    
   case {'cmpdur'}
    [Y1, void, lbls] = synthLoadData('cmp');
    [Y2, void, lbls] = synthLoadData('dur');
    Y = [Y1, Y2];

   case {'grid_vowels', ...         
         'grid_consonants', ...
         'grid_approximants', ...
         'grid_uplosives', ...
         'grid_vplosives', ...
         'grid_plosives', ...
         'grid_nasals', ...
         'grid_ufricatives', ...
         'grid_vfricatives', ...
         'grid_fricatives', ...
         'grid_affricatives', ...
         'grid_silence', ...
         'grid_other', ...
         'grid_all'}
    try 
      load([baseDir dataset '.mat']);
    catch
      
      fullList = {'iy2', 'iy4', 'iy5', 'iy6', ...     % Vowels (iy, ih, eh, ey, ae, aa, aw, ah, ao, ow, uw, ax, ia)
                'ih1', 'ih3', 'ih5', ...
                'eh1', 'eh2', 'eh4', 'eh5', 'eh6', ...
                'ey1', 'ey4', 'ey5', ...
                'ae3', ...
                'aa4', ...
                'aw6', ...
                'ay2', 'ay3', 'ay4', 'ay5', ...
                'ah5', ...
                'ao5', ...
                'ow4', 'ow5', ...
                'uw2', 'uw4', 'uw5', 'uw6', ...
                'ax4', 'ax6', ...
                'ia5', ...
                'w2', 'w3', 'w4', 'w5', ... % approximants -- kind of between a vowel and consonant (w, y, r, l)
                'y4', ...
                'r2', 'r5', ...
                'l1', 'l2', 'l4', 'l6', ... 
                'b1', 'b2', 'b3', 'b4', ... % Consonants  - plosives, voiced (b, d, g)
                'd2', 'd4', ...
                'g2', 'g6', ...
                'p1', 'p4', 'p6', ... % Consonants - plosives, unvoiced (p, t, k)
                't1', 't2', 't3', 't4', 't5',...
                'k4', 'k5', ...
                'm4', ... % Consonants - nasals (m, n)
                'n1', 'n2', 'n3', 'n4', 'n5', 'n6', ...
                'f4', 'f5', ... % Consonants - fricatives, unvoiced (f, s, th, sh)
                's1', 's4', 's5', 's6', ...
                'th5', ...
                'v4', 'v5', ... % Consonant - fricatives, voiced (v, z, dh, zh)
                'z4', 'z5', 'z6', ...
                'dh3', ...
                ...%(but consonsant can also be grouped strong/weak:   strong = s, sh, z, zh;  weak=f,v,th,dh)
                'ch4', ...  % Consonant - affricative  (basically a blend between a plosive and a fricative) (ch  =  t + sh, jh = d + zh)
                'jh4', ...
                'sil', ... %silence
                'sp'}; % ask Jon

      vowelInd = 1:33;
      consonantInd = 34:84;
      approximantInd = 34:44;
      vplosivesInd = 45:52;
      uplosivesInd = 53:62;
      nasalsInd = 63:69;
      ufricativesInd = 70:76;
      vfricativesInd = 77:82;
      affricativesInd = 83:84;
      silenceInd = 85;
      otherInd = 86;
      constanants = {};
      all = {};
      switch dataset(6:end)
       case 'vowels'
        phoneList = fullList(vowelInd);
       case 'consonants'
        phoneList = fullList(consonantInd);
       case 'approximants'
        phoneList = fullList(approximantInd);
       case 'uplosives'
        phoneList = fullList(uplosivesInd);
       case 'vplosives'
        phoneList = fullList(vplosivesInd);
       case 'plosives'
        phoneList = fullList([uplosivesInd vplosivesInd]);
       case 'nasals'
        phoneList = fullList(nasalsInd);
       case 'ufricatives'
        phoneList = fullList(ufricativesInd);
       case 'vfricatives'
        phoneList = fullList(vfricativesInd);
       case 'fricatives'
        phoneList = fullList([ufricativesInd vfricativesInd]);
       case 'affricatives'
        phoneList = fullList(affricativesInd);
       case 'silence'
        phoneList = fullList(silenceInd);
       case 'other'
        phoneList = fullList(otherInd);
       case 'all'
        phoneList = fullList;
      end
      
      motherDir = '/local/data/datasets/synthesis/models/';
      fileLoc = '/ver3/cmp/monophoneC.mmf';
      numId = 34;
      i = 1;
      idStr = ['id' num2str(i)];
      fileName = [motherDir idStr fileLoc];
      [meanValCell, varValCell] = htkLoadMmf(phoneList, fileName);
      meanVals = [meanValCell{:, :}];
      varVals = [varValCell{:, :}];
      meanVals = zeros(numId, size(meanVals, 2));
      varVals = zeros(numId, size(varVals, 2));
      
      for i = 1:numId
        idStr = ['id' num2str(i)];
        fileName = [motherDir idStr fileLoc];
        [meanValCell, varValCell] = htkLoadMmf(phoneList, fileName);
        meanVals(i, :) = [meanValCell{:, :}];
        varVals(i, :) = [varValCell{:, :}];
      end
      meanVals = meanVals./sqrt(varVals);
      varVals = 0.5*log(varVals);
      meanMean = mean(meanVals);
      meanLogStd = mean(varVals);
      for i = 1:numId 
        meanVals(i, :) = meanVals(i, :) - meanMean;
        varVals(i, :) = varVals(i, :) - meanLogStd;
      end
      varMeanVals = var(meanVals(:));
      varVarVals = var(varVals(:));
      Y = [meanVals/sqrt(varMeanVals) varVals/sqrt(varVarVals)];
      lbls = [1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, ...
              0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0]';
      lbls = [lbls ~lbls];
      save([baseDir dataset '.mat'], 'Y', 'lbls', 'meanMean', 'meanLogStd', 'varMeanVals', 'varVarVals');
    end
    if size(lbls, 2) == 1
      lbls = [lbls ~lbls];
      save([baseDir dataset '.mat'], 'Y', 'lbls', 'meanMean', 'meanLogStd', 'varMeanVals', 'varVarVals');
    end
    retLbls = cell(2,1);
    retLbls{1} = lbls;
    %/~retLbls{2} = {'JB', 'Martin Cooke', 'Stuart Cunningham', 'Sue Harding', ...
%                  'Jonny Laidler', 'Stuart Wriggley', 'Gillian', 'Mike Stannett', ...
%                  'John Karn', 'UG', 'Anna', 'UG', ...
%                  'UG', 'UG', 'Female', 'UG', ...
%                  'UG', 'UG',  'UG', 'UG', ...
%                  'UG', 'UG', 'UG', 'UG', ...
%                  'UG', 'James Carmichael', 'Rob Mill', 'Matt Gibson', ...
%                  'Sarah Simpson', 'Vinny Wan', 'Lucy', 'NDL', ...
%                  'UG', 'UG'};
                  %~/
                  retLbls{2} = {'JB', 'MPC', 'Stuart Cunningham', 'Sue Harding', ...
                                'Jonny Laidler', 'Stuart Wriggley', 'South Yorkshire G', 'MS', ...
                                'JK', 'UG 1', 'South Yorkshire A', 'UG 2', ...
                                'UG 3', 'UG 4', 'UG 5', 'UG 6', ...
                                'UG 7', 'UG 8',  'UG 9', 'UG 10', ...
                                'UG 11', 'UG 12', 'UG 13', 'UG14', ...
                                'UG', 'West Indies', 'Rob Mill', 'Scottish', ...
                                'Sarah Simpson', 'Vinny Wan', 'South Yorkshire L', 'NDL', ...
                                'UG 15', 'UG 16'};
    lbls = retLbls;

   case 'cmu35WalkJog'
    try 
      load([baseDir 'cmu35WalkJog.mat']);
    catch
      [void, errid] = lasterr;
      if strcmp(errid, 'MATLAB:load:couldNotReadFile');
        skel = acclaimReadSkel([baseDir 'mocap' dirSep 'cmu' dirSep '35' dirSep '35.asf']);
        examples = ...
            {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', ...
             '11', '12', '13', '14', '15', '16', '17', '19', '20', ...
             '21', '22', '23', '24', '25', '26', '28', '30', '31', '32', '33', '34'};
        testExamples = {'18', '29'};
        % Label differently for each sequence
        exlbls = eye(31);
        testexlbls = eye(2);
        totLength = 0;
        totTestLength = 0;
        for i = 1:length(examples)
          [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep '35' dirSep '35_' ...
                              examples{i} '.amc'], skel);
          tY{i} = tmpchan(1:4:end, :);
          tlbls{i} = repmat(exlbls(i, :), size(tY{i}, 1), 1);
          totLength = totLength + size(tY{i}, 1);
        end
        Y = zeros(totLength, size(tY{1}, 2));
        lbls = zeros(totLength, size(tlbls{1}, 2));
        endInd = 0;
        for i = 1:length(tY)
          startInd = endInd + 1;
          endInd = endInd + size(tY{i}, 1);
          Y(startInd:endInd, :) = tY{i};
          lbls(startInd:endInd, :) = tlbls{i};
        end
        for i = 1:length(testExamples)
          [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep '35' dirSep '35_' ...
                              testExamples{i} '.amc'], skel);
          tYtest{i} = tmpchan(1:4:end, :);
          tlblstest{i} = repmat(testexlbls(i, :), size(tYtest{i}, 1), 1);
          totTestLength = totTestLength + size(tYtest{i}, 1);
        end
        Ytest = zeros(totTestLength, size(tYtest{1}, 2));
        lblstest = zeros(totTestLength, size(tlblstest{1}, 2));
        endInd = 0;
        for i = 1:length(tYtest)
          startInd = endInd + 1;
          endInd = endInd + size(tYtest{i}, 1);
          Ytest(startInd:endInd, :) = tYtest{i};
          lblstest(startInd:endInd, :) = tlblstest{i};
        end
        save([baseDir 'cmu35WalkJog.mat'], 'Y', 'lbls', 'Ytest', 'lblstest');
      else
        error(lasterr);
      end
    end

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
    lbls = 'connect';
    
   case 'brendan'
    load([baseDir 'frey_rawface.mat']);
    Y = double(ff)';

    
   case 'usps'
    
    % Fix seeds
    randn('seed', 1e5);
    rand('seed', 1e5);

    load([baseDir 'usps_train.mat']);
    [ALL_T, sortIndices] = sort(ALL_T);
    ALL_DATA = ALL_DATA(sortIndices(:), :);
    Y = ALL_DATA;
    lbls = zeros(size(ALL_T,1),10);
    for digit = 0:9;
      ind1 = min(find(ALL_T==digit));
      ind2 = max(find(ALL_T==digit));
      lbls(ind1:ind2,digit+1) = ones(size(ind1:ind2,2),1);
    end
    
    % Fix seeds
    randn('seed', 1e5);
    rand('seed', 1e5);
    
    load([baseDir 'usps_test.mat']);
    [ALL_T, sortIndices] = sort(ALL_T);
    ALL_DATA = ALL_DATA(sortIndices(:), :);
    Ytest = ALL_DATA;
    lblstest = zeros(size(ALL_T,1),10);
    for digit = 0:9;
      ind1 = min(find(ALL_T==digit));
      ind2 = max(find(ALL_T==digit));
      lblstest(ind1:ind2,digit+1) = ones(size(ind1:ind2,2),1);
    end
    
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

   case 'swissRollFull'
    load([baseDir 'swiss_roll_data']);
    Y = X_data(:, 1:3000)';

    %/~
   case 'fourWalks'
    try 
      load([baseDir 'fourWalks.mat']);
    catch
      [void, errid] = lasterr;
      if strcmp(errid, 'MATLAB:load:couldNotReadFile');
        
        subj = {'35', '10', '12', '16'};
        examples = {'02', '04', '01', '15'};
        exlbls = eye(4);
        startVals = [55, 222, 22, 62];
        endVals = [338, 499, 328, 342];
        diffInd = [1 3 5];
        totLength = 0;
        totTestLength = 0;
        for i = 1:length(examples)
          skel = acclaimReadSkel([baseDir 'mocap' dirSep 'cmu' dirSep subj{i} dirSep subj{i} ...
                              '.asf']);
          [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep subj{i} dirSep subj{i} '_' ...
                              examples{i} '.amc'], skel);
          tY{i} = tmpchan(startVals(i):4:endVals(i), :);
          % Make given indices differences.
          tY{i}(:, diffInd) = tmpchan(startVals(i)+1:4:endVals(i)+1, diffInd)...
              - tY{i}(:, diffInd);
          tlbls{i} = repmat(exlbls(i, :), size(tY{i}, 1), 1);
          totLength = totLength + size(tY{i}, 1);
        end
        Y = zeros(totLength, size(tY{1}, 2));
        lbls = zeros(totLength, size(tlbls{1}, 2));
        endInd = 0;
        for i = 1:length(tY)
          startInd = endInd + 1;
          endInd = endInd + size(tY{i}, 1);
          Y(startInd:endInd, :) = tY{i};
          lbls(startInd:endInd, :) = tlbls{i};
        end
        Ytest = [];
        lblstest = [];
        subj = {'35', '12', '16', '12', '07', '07', '08', '08'};
        examples = {'03', '02', '21', '03', '01', '02', '01', '02'};
        testexlbls = eye(8);
        startVals = ones(1, 8);
        endVals = 200*ones(1, 8);
        totLength = 0;
        totTestLength = 0;
        for i = 1:length(examples)
          skel = acclaimReadSkel([baseDir 'mocap' dirSep 'cmu' dirSep subj{i} dirSep subj{i} ...
                              '.asf']);
          [tmpchan, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep subj{i} dirSep subj{i} '_' ...
                              examples{i} '.amc'], skel);
          tYtest{i} = tmpchan(startVals(i):4:endVals(i), :);
          tYtest{i}(:, diffInd) = tmpchan(startVals(i)+1:4:endVals(i)+1, diffInd)...
              - tYtest{i}(:, diffInd);
          tlblstest{i} = repmat(testexlbls(i, :), size(tYtest{i}, 1), 1);
          totLength = totLength + size(tYtest{i}, 1);
        end
        Ytest = zeros(totTestLength, size(tYtest{1}, 2));
        lblstest = zeros(totTestLength, size(tlblstest{1}, 2));
        endInd = 0;
        for i = 1:length(tYtest)
          startInd = endInd + 1;
          endInd = endInd + size(tYtest{i}, 1);
          Ytest(startInd:endInd, :) = tYtest{i};
          lblstest(startInd:endInd, :) = tlblstest{i};
        end
        save([baseDir 'fourWalks.mat'], 'Y', 'lbls', 'Ytest', 'lblstest');
      else
        error(lasterr);
      end
    end

   case 'walkSitJog'
    Y = mocapLoadTextData([baseDir 'walkSitJog']);
    Y = Y(1:4:end, :)*200;
    

   case 'robotAtrium'
    Y = parseWirelessData([baseDir 'uw-atrium.txt']);
    Y = (Y + 85)/15;
    Y = Y(1:281, :);
    lbls = 'connect'
   
    Ytest = parseWirelessData([baseDir 'uw-atrium.txt']);
    Ytest = (Ytest + 85)/15;
    Ytest = Ytest(282:end, :);

    
   case 'accWalkSitJog'
    try
      load([baseDir 'walkSitJog.mat']); 
    catch    
      [void, errid] = lasterr;
      if strcmp(errid, 'MATLAB:load:couldNotReadFile')
        skel = acclaimReadSkel([baseDir 'mocap' dirSep 'cmu' dirSep '86' dirSep '86.asf']);
        [channels, skel] = acclaimLoadChannels([baseDir 'mocap' dirSep 'cmu' dirSep '86' dirSep '86_10.amc'], skel);
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
    
    
   otherwise
    error('Unknown data set requested.')
    
  end
end
