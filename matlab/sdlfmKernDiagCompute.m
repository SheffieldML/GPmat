function k = sdlfmKernDiagCompute(sdlfmKern, t, covIC, type)

% SDLFMKERNDIAGCOMPUTE Compute diagonal of a SDLFM kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the switching
% dynamical latent force model kernel given a design matrix of inputs.
% ARG sdlfmKern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% ARG type : specifies the type of diagonal kernel to compute. Options are
% 'PosPos' (default), 'VelVel' and 'AccelAccel'.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : sdlfmKernParamInit, kernDiagCompute, kernCreate
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    type = 'PosPos';
end

switch type
    case 'PosPos'
        fhandle = 'sdlfmKernDiagComputeBlock';
    case 'VelVel'
        fhandle = 'sdlfmvKernDiagComputeBlock';
    case 'AccelAccel'
        fhandle = 'sdlfmaKernDiagComputeBlock';    
end

fhandle = str2func(fhandle);


%Form the basic kernels
lfmKern = struct();

% Create structures that will make easy the computation of the kernels
%spVector = [sdlfmKern.switchingTimes t(end)+0.1]; % For the last interval include the last point

spVector = [cumsum(sdlfmKern.switchingTimes) t(end)+50];

dim = zeros(1, sdlfmKern.nIntervals);

for i=1:sdlfmKern.nIntervals
    for j =1:sdlfmKern.nlfPerInt
        % Create the appropriate set of kernel structures
        lfmKern(i,j).mass = sdlfmKern.mass;
        lfmKern(i,j).spring = sdlfmKern.spring;
        lfmKern(i,j).damper = sdlfmKern.damper;
        lfmKern(i,j).inverseWidth = sdlfmKern.inverseWidth(j,i);
        lfmKern(i,j).sensitivity = sdlfmKern.sensitivity(j,i);
        lfmKern(i,j).limit = spVector(i+1) - spVector(i);
        lfmKern(i,j).isNormalised = sdlfmKern.isNormalised;
    end
    newt = t(t> spVector(i) & t<spVector(i+1));
    dim(i) = length(newt);
end


kyy = zeros(sdlfmKern.nIntervals,1);kyv = zeros(sdlfmKern.nIntervals,1);
kvy = zeros(sdlfmKern.nIntervals,1);kvv = zeros(sdlfmKern.nIntervals,1);
kyy(1,1) = covIC(1,1); 
kvy(1,1) = covIC(2,1); 
kyv(1,1) = covIC(1,2);
kvv(1,1) = covIC(2,2);

k = zeros(sum(dim),1);

startValOne = 1;
endValOne   = 0;

for i=1:sdlfmKern.nIntervals
    endValOne = endValOne + dim(i);
    k(startValOne:endValOne) = fhandle(lfmKern(i,:), ...
        t(startValOne:endValOne) - spVector(i), kyy(i), kyv(i), kvy(i), kvv(i));
    tInit = spVector(i+1) - spVector(i);
    kyy(i+1) = sdlfmKernDiagComputeBlock(lfmKern(i,:), ...
        tInit, kyy(i), kyv(i), kvy(i), kvv(i));
    % VEL -- POS
    kvy(i+1) = sdlfmvXsdlfmKernComputeBlock(lfmKern(i,:), ...
        lfmKern(i,:), tInit, tInit , ...
        kyy(i), kyv(i), kvy(i), kvv(i), i);
    % POS -- VEL
    kyv(i+1) = sdlfmXsdlfmvKernComputeBlock(lfmKern(i,:), ...
        lfmKern(i,:), tInit, tInit, ...
        kyy(i), kyv(i), kvy(i), kvv(i), i);
    % VEL -- VEL
    kvv(i+1) = sdlfmvKernDiagComputeBlock(lfmKern(i,:), ...
        tInit, kyy(i), kyv(i), kvy(i), kvv(i));
    startValOne = endValOne + 1;
end

k = real(k);
