% DEMHAND5 Oil data with fully independent training conditional.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'hand';
experimentNo = 5;



load points_all
p_some = p_all(1:10:end);

y = p_all{1};
y(:, 1) = y(:, 1)/sqrt(var(y(:, 1)));
y(:, 2) = y(:, 2)/sqrt(var(y(:, 2)));

t = linspace(0, 360, length(y)+1)';
t = t(1:end-1);
periodicBase =[1 2 3 5 11 22];
numPeriodics = length(periodicBase);
baset = t;
basey = y;

starty = basey(end, :);
startt = baset(end, :)-360;
pts = 3;
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
t = deg2rad(t);

tnsKern{1} = 'tensor';
for i = 1:numPeriodics 
  tnsKern{i+1} = 'rbfperiodic';
end
kern = kernCreate(t, {'cmpnd', ...
                    tnsKern, 'bias', 'white'});
kernPeriod = kern.comp{1};

for i = 1:length(kernPeriod.comp)
  kernPeriod.comp{i}.period = 2*pi/periodicBase(i);
end
kern.comp{1} = kernPeriod;

options = gpOptions('ftc');
options.kern = kern;
model = gpCreate(1, 2, t, y, ...
                         options);
model.optimiser = 'conjgrad';
model = gpOptimise(model, 1, 200);
