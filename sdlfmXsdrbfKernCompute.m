function [K, Kc] = sdlfmXsdrbfKernCompute(sdlfmKern, sdrbfKern, t1, t2, type)

% SDLFMXSDRBFKERNCOMPUTE Cross kernel between a SDLFM and a SDRBF kernels.
% FORMAT
% DESC computes cross kernel terms between a switching dynamical
% LFM kernels and a switching dynamical RBF kernel for the multiple output 
% kernel.
% ARG sdlfmKern : the kernel structure associated with the SDLFM
% kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix. K is a cell where the
% number of cells is related to the number of latent forces per interval
% RETURN Kc : the kernel matrices grouped in cells. Suitable form for the
% sparse approximations
%
% FORMAT
% DESC computes cross kernel terms between a SDLFM and a SDRBF kernels for
% the multiple output kernel.
% ARG sdlfmKern : the kernel structure associated with the SDLFM
% kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix. K is a cell where the
% number of cells is related to the number of latent forces per interval
% RETURN Kc : the kernel matrices grouped in cells. Suitable form for the
% sparse approximations
%
% FORMAT
% DESC computes cross kernel terms a SDLFM and a SDRBF kernels for
% the multiple output kernel. The SDLFM kernel can correspond to Position
% (default), Velocity or Acceleration. The type of kernel to be computed 
% is specified in 'type'. 
% ARG sdlfmKern : the kernel structure associated with the SDLFM
% kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG type : specifies the type of kerne to be computed
% RETURN K : block of values from kernel matrix. K is a cell where the
% number of cells is related to the number of latent forces per interval
% RETURN Kc : the kernel matrices grouped in cells. Suitable form for the
% sparse approximations
%
% SEEALSO : sdlfmKernParamInit, sdlfmKernCompute, sdlfmKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 5
    type = 'Pos';
    if nargin < 4
        t2 = t1;
    end
end

if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end

% compareInverseWidth = sdlfmKern.inverseWidth == sdrbfKern.inverseWidth;
% if sum(sum(compareInverseWidth))~=(sdlfmKern.nIntervals*sdlfmKern.nlfPerInt)
%     error('Kernels cannot be cross combined if they have different inverse widths.')
% end
% compareSwitchingTimes = sdlfmKern.switchingTimes == sdrbfKern.switchingTimes;
% if sum(sum(compareSwitchingTimes))~=sdlfmKern.nIntervals
%     error('Kernels cannot be cross combined if they have different switching points.')
% end

switch type
    case 'Pos'
        fhandle = 'sdlfmXsdrbfKernComputeBlock';
    case 'Vel'
        fhandle = 'sdlfmvXsdrbfKernComputeBlock';
    case 'Accel'
        fhandle = 'sdlfmaXsdrbfKernComputeBlock';   
end

fhandle = str2func(fhandle);
    
%Form the basic kernels
lfmKern = struct();
rbfKern = struct();

% Create structures that will make easy the computation of the kernels
%spVector = [sdlfmKern.switchingTimes t1(end)+0.1]; % For the last interval include the last point

spVector = [cumsum(sdlfmKern.switchingTimes) t1(end)+50];

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

if sum(dim1)~=length(t1) || sum(dim2)~=length(t2)
    error('A problem with the dimensions of the switching intervals occured')
end

% Compute some necessary constants

generalConst = sdlfmKernComputeConstant(sdlfmKern.nIntervals, ...
    lfmKern(1,1), lfmKern(1,1), spVector);

K = cell(sdrbfKern.nlfPerInt,1);
Kc = cell(sdrbfKern.nlfPerInt,1);


for q=1:sdrbfKern.nlfPerInt
    startValOne = 1;
    endValOne   = 0;
    tK = zeros(sum(dim1), sum(dim2));
    for i=1:sdlfmKern.nIntervals
        endValOne = endValOne + dim1(i);
        startValThree = 1;
        endValThree = 0;
        for j=1:i
            if i>j
                lfmKernLocal = lfmKern(j,q);
                rbfKernLocal = rbfKern(j,q);
            else
                lfmKernLocal = lfmKern(i,q);
                rbfKernLocal = rbfKern(i,q);
            end
            endValThree = endValThree + dim2(j);
            tK(startValOne:endValOne, startValThree:endValThree) = fhandle(lfmKernLocal, ...
                rbfKernLocal, t1(startValOne:endValOne) - spVector(i), ...
                t2(startValThree:endValThree) - spVector(j), i, j, generalConst);
            startValThree = endValThree + 1;
        end
        startValOne = endValOne + 1;
    end
    tK = real(tK);
    K{q} = tK;
end

