% DEMHAND1 Oil data with fully independent training conditional.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'hand';
experimentNo = 1;


numPoints = 100;
numBumps = 15;
t = [155 215 290];
     
littleFinger =  linspace(315, 335, 8);
littleFinger = littleFinger(1:end-1);
ringFinger = linspace(335, 360, 8);
ringFinger = ringFinger(1:end-1);
middleFinger = linspace(0, 25, 8);
middleFinger = middleFinger(1:end-1);
indexFinger = linspace(25, 50, 8);
thumb = linspace(75, 105, 6);

t  = [t littleFinger ringFinger middleFinger indexFinger thumb];
t = deg2rad(t)';


numPeriodics = 21;
tnsKern{1} = 'tensor';
for i = 1:numPeriodics 
  tnsKern{i+1} = 'rbfperiodic';
end
kern = kernCreate(t, {'cmpnd', ...
                    tnsKern, 'bias', 'white'});
kernPeriod = kern.comp{1};

for i = 1:length(kernPeriod.comp)
  kernPeriod.comp{i}.period = 2*pi/i;
end
kern.comp{1} = kernPeriod;

sigma = [5 5 10 5 5]/5;
sigma = sigma/3;
centers = [10 32 90 330 350];
centers = [centers+sigma centers centers-sigma];
sigma = repmat(sigma, 1, 3);


%kernPeriod.nParams = kernPeriod.lengthScaleFunc.numParams + 1;
%kern.comp{1} = kernPeriod;
%kern.nParams = 0;
%for i = 1:length(kern.comp)
%  kern.nParams = kern.nParams + kern.comp{i}.nParams;
%end
%kern.paramGroups = speye(kern.nParams);
load points_all
p_some = p_all(1:10:end);

baset = deg2rad(t);
basey = p_all{1};
starty = basey(end, :);
startt = baset(end, :);
pts = 5;
spacing = linspace(0, 1, pts + 1)';
spacing = spacing(1:end-1);
y = [];
t = [];
for i = 1:size(basey, 1)
  yvec = basey(i, :) - starty;
  tvec = baset(i) - startt;
  y = [y; spacing*yvec+repmat(starty, pts, 1)];
  t = [t; spacing*tvec+repmat(startt, pts, 1)];
  starty = basey(i, :);
  startt = baset(i, :);
end

options = gpOptions('ftc');
options.kern = kern;
model = gpCreate(1, 2, t, y, ...
                         options);
model.optimiser = 'conjgrad';
model = gpOptimise(model, 1, 2000);
