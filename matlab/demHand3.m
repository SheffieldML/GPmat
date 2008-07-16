load points_all
p_some = p_all(1:10:end);

dataSetName = 'hand';
experimentNo = 3;


numPoints = 100;
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
[t, ind] = sort(t);

baset = deg2rad(t)';
load points_all
indY = 1:150:length(p_all);
basey = [];
for i = 1:length(indY);
  basey = [basey p_all{indY(i)}(ind, :)];
end
starty = basey(end, :);
startt = baset(end, :)-2*pi;
pts = 8;
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





kern = kernCreate(t, {'cmpnd', 'gibbsperiodic', 'bias', 'white'});
kernGibbs = kern.comp{1};


sigma = [5 5 10 5 5];
sigma = sigma/3;
centers = [10 32 90 330 350];
centers = [centers+sigma centers centers-sigma];
sigma = repmat(sigma, 1, 3);
numBumps = length(centers);

options = rbfperiodicOptions(numBumps);
kernGibbs.lengthScaleFunc = modelCreate('rbfperiodic', kernGibbs.inputDimension, 1, options);

kernGibbs.lengthScaleFunc.sigma2 = deg2rad(sigma).*deg2rad(sigma);
kernGibbs.lengthScaleFunc.thetaBar = deg2rad(centers);
kernGibbs.lengthScaleFunc.bias = -2;
kernGibbs.lengthScaleFunc.weights = -1*ones(numBumps, 1);

kernGibbs.nParams = kernGibbs.lengthScaleFunc.numParams + 1;
kern.comp{1} = kernGibbs;
kern.comp{1}.variance = 1071.071;
kern.comp{2}.variance = 286;
kern.comp{3}.variance = 0.1491;
kern.nParams = 0;
for i = 1:length(kern.comp)
  kern.nParams = kern.nParams + kern.comp{i}.nParams;
end
kern.paramGroups = speye(kern.nParams);


options = gpOptions('ftc');
options.kern = kern;
for i = 1:size(y, 2);
  y(:, i) = y(:, i)/sqrt(var(y(:, i)));
end
model = gpCreate(1, size(y, 2), t, y, ...
                         options);
model = gpOptimise(model, 1, 100);
