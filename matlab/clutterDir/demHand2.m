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


kern = kernCreate(t, {'cmpnd', ...
                    {'tensor', 'rbfperiodic', 'rbfperiodic', 'rbfperiodic', ...
                    'rbfperiodic', 'rbfperiodic', 'rbfperiodic', 'rbfperiodic'}, ...
                    'white'});
kernPeriod = kern.comp{1};

kernPeriod.comp{2}.period = pi;
kernPeriod.comp{2}.period = pi/2;
kernPeriod.comp{3}.period = pi/3;
kernPeriod.comp{4}.period = pi/4;
kernPeriod.comp{5}.period = pi/5;
kernPeriod.comp{6}.period = pi/6;
kernPeriod.comp{7}.period = pi/7;
kernPeriod.comp{8}.period = pi/7;
kernPeriod.comp{9}.period = pi/7;
kernPeriod.comp{10}.period = pi/7;
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


Y = p_some;

for i = 1:length(Y)
  X{i} = t;
  Y{i} = (Y{i}-140)/33;
end
options = multimodelOptions('gp', length(Y), 'ftc');
options.compOptions.kern = kern;
model = multimodelCreate(1, 2, X, Y, ...
                         options);

%model = modelOptimise(model, 1, 10);

% Save the results.
%capName = dataSetName;;
%capName(1) = upper(capName(1));
%save(['dem' capName num2str(experimentNo) '.mat'], 'model');
