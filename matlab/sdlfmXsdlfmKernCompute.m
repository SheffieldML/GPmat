function K = sdlfmXsdlfmKernCompute(sdlfmKern1, sdlfmKern2, t1, t2, covIC, type)

% SDLFMXSDLFMKERNCOMPUTE Compute a cross kernel between two SDLFM kernels.
% FORMAT
% DESC computes cross kernel terms between two switching dynamical
% LFM kernels for the multiple output kernel.
% ARG sdlfmKern1 : the kernel structure associated with the first SDLFM
% kernel.
% ARG sdlfmKern2 : the kernel structure associated with the second SDLFM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% ARG covIC : covariance for the initial conditions
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two SDLFM kernels for
% the multiple output kernel.
% ARG sdlfmKern1 : the kernel structure associated with the first SDLFM
% kernel.
% ARG sdlfmKern2 : the kernel structure associated with the second SDLFM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covIC : covariance for the initial conditions
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two SDLFM kernels for
% the multiple output kernel. The SDLFM kernels can correspond to Position
% X Position (default), Velocity X Position, Velocity X Velocity,
% Acceleration X Position, Acceleration X Velocity, Acceleration X
% Acceleration. The type of kernel to be computed is specified in 'type'. 
% ARG sdlfmKern1 : the kernel structure associated with the first SDLFM
% kernel.
% ARG sdlfmKern2 : the kernel structure associated with the second SDLFM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covIC : covariance for the initial conditions
% ARG type : specifies the type of kerne to be computed
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : sdlfmKernParamInit, sdlfmKernCompute, sdlfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 6
    type = 'PosPos';
    if nargin < 5
        covIC = t2;
        t2 = t1;
    end
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
        fhandle = 'sdlfmXsdlfmKernComputeBlock';
    case 'VelPos'
        fhandle = 'sdlfmvXsdlfmKernComputeBlock';
    case 'VelVel'
        fhandle = 'sdlfmvXsdlfmvKernComputeBlock';
    case 'AccelPos'
        fhandle = 'sdlfmaXsdlfmKernComputeBlock';
    case 'AccelVel'
        fhandle = 'sdlfmaXsdlfmvKernComputeBlock';
    case 'AccelAccel'
        fhandle = 'sdlfmaXsdlfmaKernComputeBlock';    
end

fhandle = str2func(fhandle);
    
%Form the basic kernels
lfmKern1 = struct();
lfmKern2 = struct();

% Create structures that will make easy the computation of the kernels

spVector = [cumsum(sdlfmKern1.switchingTimes) t1(end)+50];

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

if sum(dim1)~=length(t1) || sum(dim2)~=length(t2)
    error('A problem with the dimensions of the switching intervals occured')
end

% Compute some necessary constants

generalConst = sdlfmKernComputeConstant(sdlfmKern1.nIntervals, ...
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

startValOne = 1;
endValOne   = 0;

for i=1:sdlfmKern1.nIntervals
    endValOne = endValOne + dim1(i);
    startValThree = 1;
    endValThree = 0;
    for j=1:sdlfmKern1.nIntervals
        if i>j
            lfmKern1Local = lfmKern1(j,:);
            lfmKern2Local = lfmKern2(j,:);
        else
            lfmKern1Local = lfmKern1(i,:);
            lfmKern2Local = lfmKern2(i,:);
        end            
        endValThree = endValThree + dim2(j);
        % POS -- POS (Kernel and initial positions)
        K(startValOne:endValOne, startValThree:endValThree) = fhandle(lfmKern1Local, ...
            lfmKern2Local, t1(startValOne:endValOne) - spVector(i), t2(startValThree:endValThree) - spVector(j), ...
            kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, generalConst);       
        % Time vector initial conditions
        tInit1 = [spVector(i) - spVector(i);spVector(i+1) - spVector(i)];
        tInit2 = [spVector(j) - spVector(j);spVector(j+1) - spVector(j)];        
        tempPosPos{i,j} = sdlfmXsdlfmKernComputeBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, ...
            kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, generalConst);
        kyy = organizeIC(kyy, tempPosPos, i, j);
        % VEL -- POS
        tempVelPos{i,j} = sdlfmvXsdlfmKernComputeBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2 , ...
            kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, generalConst);
        kvy = organizeIC(kvy, tempVelPos, i, j);
        % POS -- VEL
        tempPosVel{i,j} = sdlfmXsdlfmvKernComputeBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, ...
            kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, generalConst);
        kyv = organizeIC(kyv, tempPosVel, i, j);
        % POS -- VEL
        tempVelVel{i,j} = sdlfmvXsdlfmvKernComputeBlock(lfmKern1Local, ...
            lfmKern2Local, tInit1, tInit2, ...
            kyy(i,j), kyv(i,j), kvy(i,j), kvv(i,j), i, j, generalConst);
        kvv = organizeIC(kvv, tempVelVel, i, j);
        startValThree = endValThree + 1;
    end
    startValOne = endValOne + 1;
end

K = real(K);








