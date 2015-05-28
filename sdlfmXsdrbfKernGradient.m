function [g1, g2] = sdlfmXsdrbfKernGradient(sdlfmKern, sdrbfKern, t1, t2, covGrad, type)

% SDLFMXSDRBFKERNGRADIENT Cross gradient between a SDLFM and a SDRBF kernels.
% FORMAT
% DESC computes the gradients of an objective function with respect
% to cross kernel terms between a switching dynamical LFM kernels and a
% switching dynamical RBF kernel for the multiple output kernel.
% ARG sdlfmKern : the kernel structure associated with the SDLFM
% kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of SDLFM kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of SDRBF kernel.
%
% FORMAT
% DESC computes the gradients of an objective function with respect
% to cross kernel terms between a switching dynamical LFM kernels and a
% switching dynamical RBF kernel for the multiple output kernel.
% ARG sdlfmKern : the kernel structure associated with the SDLFM
% kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of SDLFM kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of SDRBF kernel.
%
% FORMAT
% DESC computes the gradients of an objective function with respect
% to cross kernel terms between a switching dynamical LFM kernels and a
% switching dynamical RBF kernel for the multiple output kernel.
% The SDLFM kernel can correspond to Position (default), Velocity or
% Acceleration. The type of kernel to be computed is specified in 'type'.
% ARG sdlfmKern : the kernel structure associated with the SDLFM
% kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG type : specifies the type of kerne to be computed
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of SDLFM kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of SDRBF kernel.
%
% SEEALSO : sdlfmKernParamInit, sdlfmKernCompute, sdlfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 6
    type = 'Pos';
    if nargin < 5
        covGrad = t2;
        t2 = t1;
    end
end

if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end

compareInverseWidth = sdlfmKern.inverseWidth == sdrbfKern.inverseWidth;
if sum(sum(compareInverseWidth))~=(sdlfmKern.nIntervals*sdlfmKern.nlfPerInt)
    error('Kernels cannot be cross combined if they have different inverse widths.')
end
compareSwitchingTimes = sdlfmKern.switchingTimes == sdrbfKern.switchingTimes;
if sum(sum(compareSwitchingTimes))~=sdlfmKern.nIntervals
    error('Kernels cannot be cross combined if they have different switching points.')
end

switch type
    case 'Pos'
        fhandle = 'sdlfmXsdrbfKernGradientBlock';
    case 'Vel'
        fhandle = 'sdlfmvXsdrbfKernGradientBlock';
    case 'Accel'
        fhandle = 'sdlfmaXsdrbfKernGradientBlock';
end

fhandle1 = str2func(fhandle);

%Form the basic kernels
lfmKern = struct();
rbfKern = struct();

% Create structures that will make easy the computation of the kernels

spVector = [cumsum(sdlfmKern.switchingTimes) t1(end)+100];


dim1 = zeros(1, sdlfmKern.nIntervals); dim2 = zeros(1, sdrbfKern.nIntervals);

for i =1:sdlfmKern.nIntervals
    for j=1:sdlfmKern.nlfPerInt
        % Create the appropriate set of kernel structures
        lfmKern(i,j).mass = sdlfmKern.mass;
        lfmKern(i,j).spring = sdlfmKern.spring;
        lfmKern(i,j).damper = sdlfmKern.damper;
        lfmKern(i,j).inverseWidth = sdlfmKern.inverseWidth(j,i);
        lfmKern(i,j).sensitivity = sdlfmKern.sensitivity(j,i);
        rbfKern(i,j).variance = sdrbfKern.variance;
        rbfKern(i,j).inverseWidth = sdrbfKern.inverseWidth(j,i);
        lfmKern(i,j).limit = spVector(i+1) - spVector(i);
        rbfKern(i,j).limit = spVector(i+1) - spVector(i);
        lfmKern(i,j).isNormalised = sdlfmKern.isNormalised;
        rbfKern(i,j).isNormalised = sdrbfKern.isNormalised;
    end
    newt1 = t1(t1> spVector(i) & t1<spVector(i+1));
    newt2 = t2(t2> spVector(i) & t2<spVector(i+1));
    dim1(i) = length(newt1);
    dim2(i) = length(newt2);
end

% Compute some necessary constants

[generalConstGrad, generalConst] = sdlfmKernGradientConstant(sdlfmKern.nIntervals, ...
    lfmKern(1,1), lfmKern(1,1), spVector);

% Initialization of the vector of gradients
gradIW = zeros(sdlfmKern.nlfPerInt, sdlfmKern.nIntervals);
gradS = zeros(sdlfmKern.nlfPerInt, sdlfmKern.nIntervals);
gradSP = zeros(1, sdlfmKern.nIntervals);
g1 = zeros(1, sdlfmKern.nParams);
g2 = zeros(1, sdrbfKern.nParams);
covGrad2 = covGrad;
for q=1:sdrbfKern.nlfPerInt
    startValOne = 1;
    endValOne   = 0;
    g1P = zeros(1,3);
    gradSPTemp = zeros(1, sdlfmKern.nIntervals);
    if iscell(covGrad2)
        covGrad = covGrad2{q};
    end
    for i=1:sdlfmKern.nIntervals
        endValOne = endValOne + dim1(i);
        startValThree = 1;
        endValThree = 0;
        for j=1:i
            if i>j
                lowest = j;
                lfmKernLocal = lfmKern(j,q);
                rbfKernLocal = rbfKern(j,q);
            else
                lowest = i;
                lfmKernLocal = lfmKern(i,q);
                rbfKernLocal = rbfKern(i,q);
            end
            endValThree = endValThree + dim2(j);
            [g1Local, g2Local, g3Local] = fhandle1(lfmKernLocal, rbfKernLocal, ...
                t1(startValOne:endValOne) - spVector(i), t2(startValThree:endValThree) - spVector(j), ...
                i, j, generalConst, generalConstGrad, ...
                covGrad(startValOne:endValOne, startValThree:endValThree));
            g1P = g1P + g1Local(1:3);
            gradIW(q, lowest) = gradIW(q, lowest) + g1Local(4) + g2Local(2);
            gradS(q, lowest) = gradS(q, lowest) + g1Local(5);
            gradSPTemp(1:length(g3Local)) = gradSPTemp(1:length(g3Local)) + g3Local;
            startValThree = endValThree + 1;
        end
        startValOne = endValOne + 1;
    end
    g1(1:3) = g1(1:3) + g1P;
    gradSP = gradSP + gradSPTemp;
end

% For g2 we leave the gradients in zero

g1(sdlfmKern.inverseWidthIndx) = gradIW(:)';
tempGradSP = fliplr(gradSP);
tempGradSP = cumsum(tempGradSP);
g1(sdlfmKern.switchingTimesIndx) = fliplr(tempGradSP);
g1(sdlfmKern.sensitivityIndx) = gradS(:)';



