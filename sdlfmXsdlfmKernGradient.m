function [g1, g2, covGradLocal] = sdlfmXsdlfmKernGradient(sdlfmKern1, ...
    sdlfmKern2, t1, t2, covGrad, covIC, type)

% SDLFMXSDLFMKERNGRADIENT Gradients of cross kernel between 2 SDLFM kernels.
% FORMAT
% DESC computes a cross gradient for a cross kernel between two switching 
% dynamical LFM kernels for the multiple output kernel.
% ARG sdlfmKern1 : the kernel structure associated with the first SDLFM
% kernel.
% ARG sdlfmKern2 : the kernel structure associated with the second SDLFM
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% ARG covIC : covariance for the initial conditions
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
% RETURN covGrad : partial covariances for kyy(t_0, t_0), kyv(t_0,t_0),
% kvy(t_0,t_0) and kvv(t_0,t_0). 
%
% FORMAT
% DESC computes a cross gradient for a cross kernel between two switching 
% dynamical LFM kernels for the multiple output kernel.
% ARG sdlfmKern1 : the kernel structure associated with the first SDLFM
% kernel.
% ARG sdlfmKern2 : the kernel structure associated with the second SDLFM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% ARG covIC : covariance for the initial conditions
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
% RETURN covGrad : partial covariances for kyy(t_0, t_0), kyv(t_0,t_0),
% kvy(t_0,t_0) and kvv(t_0,t_0). 
%
% FORMAT
% DESC computes a cross gradient for a cross kernel between two switching 
% dynamical LFM kernels for the multiple output kernel. SDLFM kernels can 
% correspond to Position X Position (default), Velocity X Position, 
% Velocity X Velocity, Acceleration X Position, Acceleration X Velocity, 
% Acceleration X Acceleration. The type of kernel for which the gradients
% are obtained is specified in 'type'. 
% ARG sdlfmKern1 : the kernel structure associated with the first SDLFM
% kernel.
% ARG sdlfmKern2 : the kernel structure associated with the second SDLFM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% ARG covIC : covariance for the initial conditions
% ARG type : specifies the type of kernel to be computed
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
% RETURN covGrad : partial covariances for kyy(t_0, t_0), kyv(t_0,t_0),
% kvy(t_0,t_0) and kvv(t_0,t_0). 
%
% SEEALSO : sdlfmKernParamInit, sdlfmKernCompute, sdlfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin == 5
    covIC = covGrad;
    covGrad = t2;
    t2 = t1;
    type = 'PosPos';
elseif nargin == 6
    type = 'PosPos';
end

if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end

compareInverseWidth = sdlfmKern1.inverseWidth == sdlfmKern2.inverseWidth;
if sum(sum(compareInverseWidth))~=(sdlfmKern1.nIntervals*sdlfmKern1.nlfPerInt)
    error('Kernels cannot be cross combined if they have different inverse widths.')
end
compareSwitchingTimes = sdlfmKern1.switchingTimes == sdlfmKern2.switchingTimes;
if sum(sum(compareSwitchingTimes))~=sdlfmKern1.nIntervals
    error('Kernels cannot be cross combined if they have different switching points.')
end

switch type
    case 'PosPos'
        fhandle1 = 'sdlfmXsdlfmKernGradientBlock';
        typeIC   = {'sdlfm', 'sdlfm'};
    case 'VelPos'
        fhandle1 = 'sdlfmvXsdlfmKernGradientBlock';
        typeIC   = {'sdlfmv', 'sdlfm'};
    case 'VelVel'
        fhandle1 = 'sdlfmvXsdlfmvKernGradientBlock';
        typeIC   = {'sdlfmv', 'sdlfmv'};
    case 'AccelPos'
        fhandle1 = 'sdlfmaXsdlfmKernGradientBlock';
        typeIC   = {'sdlfma', 'sdlfm'};
    case 'AccelVel'
        fhandle1 = 'sdlfmaXsdlfmvKernGradientBlock';
        typeIC   = {'sdlfma', 'sdlfmv'};
    case 'AccelAccel'
        fhandle1 = 'sdlfmaXsdlfmaKernGradientBlock';
        typeIC   = {'sdlfma', 'sdlfma'};
end

fhandle = str2func(fhandle1);

%Form the basic kernels
lfmKern1 = struct();
lfmKern2 = struct();

% Create structures that will make easy the computation of the kernels

spVector = [cumsum(sdlfmKern1.switchingTimes) t1(end)+100];

dim1 = zeros(1, sdlfmKern1.nIntervals); dim2 = zeros(1, sdlfmKern1.nIntervals);

for i=1:sdlfmKern1.nIntervals
    for j =1:sdlfmKern1.nlfPerInt
        % Create the appropriate set of kernel structures
        lfmKern1(i,j).mass = sdlfmKern1.mass;
        lfmKern1(i,j).spring = sdlfmKern1.spring;
        lfmKern1(i,j).damper = sdlfmKern1.damper;
        lfmKern1(i,j).inverseWidth = sdlfmKern1.inverseWidth(j,i);
        lfmKern1(i,j).sensitivity = sdlfmKern1.sensitivity(j,i);
        lfmKern2(i,j).mass = sdlfmKern2.mass;
        lfmKern2(i,j).spring = sdlfmKern2.spring;
        lfmKern2(i,j).damper = sdlfmKern2.damper;
        lfmKern2(i,j).inverseWidth = sdlfmKern2.inverseWidth(j,i);
        lfmKern2(i,j).sensitivity = sdlfmKern2.sensitivity(j,i);
        lfmKern1(i,j).limit = spVector(i+1) - spVector(i);
        lfmKern2(i,j).limit = spVector(i+1) - spVector(i);
        lfmKern1(i,j).isNormalised = sdlfmKern1.isNormalised;
        lfmKern2(i,j).isNormalised = sdlfmKern2.isNormalised;
    end
    newt1 = t1(t1> spVector(i) & t1<spVector(i+1));
    newt2 = t2(t2> spVector(i) & t2<spVector(i+1));
    dim1(i) = length(newt1);
    dim2(i) = length(newt2);
end

% Compute some necessary constants and their gradients

[generalConstGrad, generalConst] = sdlfmKernGradientConstant(sdlfmKern1.nIntervals, ...
    lfmKern1(1), lfmKern2(1), spVector);

kyy = zeros(sdlfmKern1.nIntervals);kyv = zeros(sdlfmKern1.nIntervals);
kvy = zeros(sdlfmKern1.nIntervals);kvv = zeros(sdlfmKern1.nIntervals);
kyy(1,1) = covIC(1,1); 
kvy(1,1) = covIC(2,1); 
kyv(1,1) = covIC(1,2);
kvv(1,1) = covIC(2,2);


%%% Compute initial conditions for intervals (1-2), (2-2) and (2-1)

tInit1 = [spVector(1) - spVector(1);spVector(2) - spVector(1)];
tInit2 = [spVector(1) - spVector(1);spVector(2) - spVector(1)];

% Pos -- Pos 
temp = sdlfmXsdlfmKernComputeBlock(lfmKern1(1), lfmKern2(1), ...
    tInit1, tInit2, kyy(1,1), kyv(1,1), kvy(1,1), kvv(1,1), 1, 1, generalConst);
kyy(2,2) = temp(2,2); kyy(1,2) = temp(1,2); kyy(2,1)  = temp(2,1);
% Vel -- Pos
temp = sdlfmvXsdlfmKernComputeBlock(lfmKern1(1), lfmKern2(1), ...
    tInit1, tInit2, kyy(1,1), kyv(1,1), kvy(1,1), kvv(1,1), 1, 1, generalConst);
kvy(2,2) = temp(2,2); kvy(1,2) = temp(1,2); kvy(2,1)  = temp(2,1);
% Pos -- Vel
temp = sdlfmXsdlfmvKernComputeBlock(lfmKern1(1), lfmKern2(1), ...
    tInit1, tInit2, kyy(1,1), kyv(1,1), kvy(1,1), kvv(1,1), 1, 1, generalConst);
kyv(2,2) = temp(2,2); kyv(1,2) = temp(1,2); kyv(2,1)  = temp(2,1);
% Vel -- Vel
temp = sdlfmvXsdlfmvKernComputeBlock(lfmKern1(1), lfmKern2(1), ...
    tInit1, tInit2, kyy(1,1), kyv(1,1), kvy(1,1), kvv(1,1), 1, 1, generalConst);
kvv(2,2) = temp(2,2); kvv(1,2) = temp(1,2); kvv(2,1)  = temp(2,1);


tempPosPos = cell(sdlfmKern1.nIntervals);
tempVelPos = cell(sdlfmKern1.nIntervals);
tempPosVel = cell(sdlfmKern1.nIntervals);
tempVelVel = cell(sdlfmKern1.nIntervals);
g1PP = cell(sdlfmKern1.nIntervals); g2PP = cell(sdlfmKern1.nIntervals);
g3PP = cell(sdlfmKern1.nIntervals); g4PP = cell(sdlfmKern1.nIntervals);
g1PV = cell(sdlfmKern1.nIntervals); g2PV = cell(sdlfmKern1.nIntervals);
g3PV = cell(sdlfmKern1.nIntervals); g4PV = cell(sdlfmKern1.nIntervals);
g1VP = cell(sdlfmKern1.nIntervals); g2VP = cell(sdlfmKern1.nIntervals);
g3VP = cell(sdlfmKern1.nIntervals); g4VP = cell(sdlfmKern1.nIntervals);
g1VV = cell(sdlfmKern1.nIntervals); g2VV = cell(sdlfmKern1.nIntervals);
g3VV = cell(sdlfmKern1.nIntervals); g4VV = cell(sdlfmKern1.nIntervals);
gkyy1 = cell(sdlfmKern1.nIntervals); gkyy2 = cell(sdlfmKern1.nIntervals);
gkyy3 = cell(sdlfmKern1.nIntervals); gkyy4 = cell(sdlfmKern1.nIntervals);
gkyv1 = cell(sdlfmKern1.nIntervals); gkyv2 = cell(sdlfmKern1.nIntervals);
gkyv3 = cell(sdlfmKern1.nIntervals); gkyv4 = cell(sdlfmKern1.nIntervals);
gkvy1 = cell(sdlfmKern1.nIntervals); gkvy2 = cell(sdlfmKern1.nIntervals);
gkvy3 = cell(sdlfmKern1.nIntervals); gkvy4 = cell(sdlfmKern1.nIntervals);
gkvv1 = cell(sdlfmKern1.nIntervals); gkvv2 = cell(sdlfmKern1.nIntervals);
gkvv3 = cell(sdlfmKern1.nIntervals); gkvv4 = cell(sdlfmKern1.nIntervals);

startValOne = 1;
endValOne   = 0;


% Initialization of the vector of gradients
g1 = zeros(1, sdlfmKern1.nParams);
g2 = zeros(1, sdlfmKern1.nParams);
gradS1 = zeros(sdlfmKern1.nlfPerInt, sdlfmKern1.nIntervals);
gradS2 = zeros(sdlfmKern1.nlfPerInt, sdlfmKern1.nIntervals);
gradIW1 = zeros(sdlfmKern1.nlfPerInt, sdlfmKern1.nIntervals);
gradIW2 = zeros(sdlfmKern1.nlfPerInt, sdlfmKern1.nIntervals);
gradSP = zeros(1, sdlfmKern1.nIntervals);


%%% Precomputations

% Needed to compute derivatives with respect to the covariances of the
% initial conditions

covGradLocal = sdlfmKernMeanCovPartial(lfmKern1(1), lfmKern2(1), ...
    t1(1:dim1(1)) - spVector(1), t2(1:dim2(1)) - spVector(1), ...
    covGrad(1:dim1(1), 1:dim2(1)), typeIC);

% Make initialitions for some derivatives
for k=1:3
    gkyy1{1,1}{k} = 0; gkyy2{1,1}{k} = 0;
    gkyv1{1,1}{k} = 0; gkyv2{1,1}{k} = 0;
    gkvy1{1,1}{k} = 0; gkvy2{1,1}{k} = 0;
    gkvv1{1,1}{k} = 0; gkvv2{1,1}{k} = 0;
end

gkyy3{1,1} = 0; gkyv3{1,1} = 0; gkvy3{1,1} = 0; gkvv3{1,1} = 0; 
gkyy4{1,1} = zeros(2); gkyv4{1,1} = zeros(2); 
gkvy4{1,1} = zeros(2); gkvv4{1,1} = zeros(2); 


for i=1:sdlfmKern1.nIntervals
    endValOne = endValOne + dim1(i);
    startValThree = 1;
    endValThree = 0;
    for j=1:sdlfmKern1.nIntervals
        if i>j
            lfmKern1Local = lfmKern1(j,:);
            lfmKern2Local = lfmKern2(j,:);
            lowest = j;
            biggest = i;
        else
            lfmKern1Local = lfmKern1(i,:);
            lfmKern2Local = lfmKern2(i,:);
            lowest = i;
            biggest = j;
        end
        endValThree = endValThree + dim2(j);
        % POS -- POS (Kernel and initial positions)
        [g1Local, g2Local, g3Local] = fhandle(lfmKern1Local, lfmKern2Local, ...
             t1(startValOne:endValOne) - spVector(i), t2(startValThree:endValThree) - spVector(j), ...
            kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, generalConst, generalConstGrad, ...
            covGrad(startValOne:endValOne, startValThree:endValThree));
        % Time vector initial conditions
        tInit1 = [spVector(i) - spVector(i);spVector(i+1) - spVector(i)];
        tInit2 = [spVector(j) - spVector(j);spVector(j+1) - spVector(j)];        
        tempPosPos{i,j} = sdlfmXsdlfmKernComputeBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, kyy(i,j), kyv(i,j), kvy(i,j), ...
            kvv(i,j), i, j, generalConst);
        kyy = organizeIC(kyy, tempPosPos, i, j);
        [g1PP{i,j}, g2PP{i,j}, g3PP{i,j}, g4PP{i,j}] = sdlfmXsdlfmKernGradientICBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, ...
            generalConst, generalConstGrad,gkyy1, gkyy2, gkyy3, gkyy4, gkyv1, ...
            gkyv2, gkyv3, gkyv4, gkvy1, gkvy2, gkvy3, gkvy4, gkvv1, gkvv2, gkvv3, gkvv4);
        % VEL -- POS
        tempVelPos{i,j} = sdlfmvXsdlfmKernComputeBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2 , kyy(i,j), kyv(i,j), kvy(i,j), ...
            kvv(i,j), i, j, generalConst);
        kvy = organizeIC(kvy, tempVelPos, i, j);
        % Gradients 
        [g1VP{i,j}, g2VP{i,j}, g3VP{i,j}, g4VP{i,j}] = sdlfmvXsdlfmKernGradientICBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, ...
            generalConst, generalConstGrad, gkyy1, gkyy2, gkyy3, gkyy4, gkyv1, ...
            gkyv2, gkyv3, gkyv4, gkvy1, gkvy2, gkvy3, gkvy4, gkvv1, gkvv2, gkvv3, gkvv4);                            
        % POS -- VEL
        tempPosVel{i,j} = sdlfmXsdlfmvKernComputeBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, kyy(i,j), kyv(i,j), kvy(i,j), ...
            kvv(i,j), i, j, generalConst);
        kyv = organizeIC(kyv, tempPosVel, i, j);
        % Gradients
        [g1PV{i,j}, g2PV{i,j}, g3PV{i,j}, g4PV{i,j}] = sdlfmXsdlfmvKernGradientICBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, ...
            generalConst, generalConstGrad, gkyy1, gkyy2, gkyy3, gkyy4, gkyv1, ...
            gkyv2, gkyv3, gkyv4, gkvy1, gkvy2, gkvy3, gkvy4, gkvv1, gkvv2, gkvv3, gkvv4);       
        % VEL -- VEL
        tempVelVel{i,j} = sdlfmvXsdlfmvKernComputeBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, kyy(i,j), kyv(i,j), kvy(i,j), ...
            kvv(i,j), i, j, generalConst);
        kvv = organizeIC(kvv, tempVelVel, i, j);
        % Gradients
        [g1VV{i,j}, g2VV{i,j}, g3VV{i,j}, g4VV{i,j}] = sdlfmvXsdlfmvKernGradientICBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, ...
            generalConst, generalConstGrad, gkyy1, gkyy2, gkyy3, gkyy4, gkyv1, ...
            gkyv2, gkyv3, gkyv4, gkvy1, gkvy2, gkvy3, gkvy4, gkvv1, gkvv2, gkvv3, gkvv4);
        % Organize derivatives
        [gkyy1, gkyy2, gkyy3, gkyy4] = organizeDerivatives(gkyy1, gkyy2, gkyy3, gkyy4, ...
            g1PP, g2PP, g3PP , g4PP, i,j);  
        [gkyv1, gkyv2, gkyv3, gkyv4] = organizeDerivatives(gkyv1, gkyv2, gkyv3, gkyv4, ...
            g1PV, g2PV, g3PV , g4PV, i,j);        
        [gkvy1, gkvy2, gkvy3, gkvy4] = organizeDerivatives(gkvy1, gkvy2, gkvy3, gkvy4, ...
            g1VP, g2VP, g3VP , g4VP, i,j);
        [gkvv1, gkvv2, gkvv3, gkvv4] = organizeDerivatives(gkvv1, gkvv2, gkvv3, gkvv4, ...
            g1VV, g2VV, g3VV , g4VV,i,j);
        % Be aware that these derivatives are from the parameters in
        % the interval before. This is important for inversewidths,
        % switching points and sensitivities. For the other parameters
        % it is not so important.        
        [g1LIC, g2LIC, g3LIC, g4LIC] = sdlfmXsdlfmKernGradientIC(lfmKern1Local(1), lfmKern2Local(1), ...
            t1(startValOne:endValOne) - spVector(i), t2(startValThree:endValThree) - spVector(j), ...
            gkyy1{i,j}, gkyy2{i,j}, gkyy3{i,j}, gkyy4{i,j}, gkyv1{i,j}, gkyv2{i,j}, gkyv3{i,j}, gkyv4{i,j}, ...
            gkvy1{i,j}, gkvy2{i,j}, gkvy3{i,j}, gkvy4{i,j}, gkvv1{i,j}, gkvv2{i,j}, gkvv3{i,j}, gkvv4{i,j}, ...
            covGrad(startValOne:endValOne, startValThree:endValThree), typeIC);
        % Organize the normal derivatives
        g1(sdlfmKern1.outputIndx) = g1(sdlfmKern1.outputIndx) + g1Local{1} + g1LIC{1};
        g2(sdlfmKern1.outputIndx) = g2(sdlfmKern1.outputIndx) + g2Local{1} + g2LIC{1};
        if i>=2 && j>=2
            localIW = reshape(g1LIC{2},sdlfmKern1.nlfPerInt, lowest-1);
            gradIW1(:,1:lowest) = gradIW1(:,1:lowest) + [localIW  g1Local{2}'];
            localIW = reshape(g2LIC{2},sdlfmKern1.nlfPerInt, lowest-1);
            gradIW2(:,1:lowest) = gradIW2(:,1:lowest) + [localIW g2Local{2}'];
            localS = reshape(g1LIC{3},sdlfmKern1.nlfPerInt, lowest-1);
            gradS1(:,1:lowest) = gradS1(:,1:lowest) + [localS g1Local{3}'];
            localS = reshape(g2LIC{3},sdlfmKern1.nlfPerInt, lowest-1);
            gradS2(:,1:lowest) = gradS2(:,1:lowest) + [localS g2Local{3}'];
        else
            gradIW1(:,lowest) = gradIW1(:,lowest) + g1Local{2}';
            gradIW2(:,lowest) = gradIW2(:,lowest) + g2Local{2}';
            gradS1(:,lowest) = gradS1(:,lowest) + g1Local{3}';
            gradS2(:,lowest) = gradS2(:,lowest) + g2Local{3}';
        end
        gradSP(1:biggest) = gradSP(1:biggest) + g3Local + g3LIC;
        covGradLocal = covGradLocal + g4LIC;
        startValThree = endValThree + 1;
    end
    startValOne = endValOne + 1;
end

g1(sdlfmKern1.inverseWidthIndx) = gradIW1(:)' + gradIW2(:)';
g2(sdlfmKern1.inverseWidthIndx) = zeros(1, sdlfmKern1.nlfPerInt*sdlfmKern1.nIntervals);
tempGradSP = fliplr(gradSP);
tempGradSP = cumsum(tempGradSP);
g1(sdlfmKern1.switchingTimesIndx) = fliplr(tempGradSP); 
g2(sdlfmKern1.switchingTimesIndx) = zeros(1, sdlfmKern1.nIntervals);
g1(sdlfmKern1.sensitivityIndx) = gradS1(:)';
g2(sdlfmKern1.sensitivityIndx) = gradS2(:)';


function [gkyy1, gkyy2, gkyy3, gkyy4] = organizeDerivatives(gkyy1, gkyy2, gkyy3, gkyy4, ...
    g1PP, g2PP, g3PP , g4PP, i,j)

if i==1 && j==1
    gkyy1{i+1, j+1} = g1PP{i,j}{2,2}; gkyy2{i+1, j+1} = g2PP{i,j}{2,2};
    gkyy3{i+1, j+1} = g3PP{i,j}{2,2}; gkyy4{i+1, j+1} = g4PP{i,j}{2,2};
    gkyy1{i,j+1} = g1PP{i,j}{1,2}; gkyy2{i,j+1} = g2PP{i,j}{1,2};
    gkyy3{i,j+1} = g3PP{i,j}{1,2}; gkyy4{i,j+1} = g4PP{i,j}{1,2};
    gkyy1{i+1,j} = g1PP{i,j}{2,1}; gkyy2{i+1,j} = g2PP{i,j}{2,1};
    gkyy3{i+1,j} = g3PP{i,j}{2,1}; gkyy4{i+1,j} = g4PP{i,j}{2,1}; 
end
if i==1 && j~=1
    gkyy1{i+1, j+1} = g1PP{i,j}{2,2}; gkyy2{i+1, j+1} = g2PP{i,j}{2,2};
    gkyy3{i+1, j+1} = g3PP{i,j}{2,2}; gkyy4{i+1, j+1} = g4PP{i,j}{2,2};
    gkyy1{i,j+1} = g1PP{i,j}{1,2}; gkyy2{i,j+1} = g2PP{i,j}{1,2};
    gkyy3{i,j+1} = g3PP{i,j}{1,2}; gkyy4{i,j+1} = g4PP{i,j}{1,2};
end
if i~=1 && j==1
    gkyy1{i+1, j+1} = g1PP{i,j}{2,2}; gkyy2{i+1, j+1} = g2PP{i,j}{2,2};
    gkyy3{i+1, j+1} = g3PP{i,j}{2,2}; gkyy4{i+1, j+1} = g4PP{i,j}{2,2};
    gkyy1{i+1,j} = g1PP{i,j}{2,1}; gkyy2{i+1,j} = g2PP{i,j}{2,1};
    gkyy3{i+1,j} = g3PP{i,j}{2,1}; gkyy4{i+1,j} = g4PP{i,j}{2,1}; 
end
if i~=1 && j~=1
    gkyy1{i+1, j+1} = g1PP{i,j}{2,2}; gkyy2{i+1, j+1} = g2PP{i,j}{2,2};
    gkyy3{i+1, j+1} = g3PP{i,j}{2,2}; gkyy4{i+1, j+1} = g4PP{i,j}{2,2};
end
