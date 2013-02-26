function K = sdlfmXsdlfmvKernComputeBlock(lfmKern1, lfmKern2, t1, t2, ...
    kyy, kyv, kvy, kvv, i, j, generalConst)

% SDLFMXSDLFMVKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
%
%	Description:
%
%	K = SDLFMXSDLFMVKERNCOMPUTEBLOCK(LFMKERN1, LFMKERN2, T1, T2, KYY,
%	KYV, KVY, KVV, I, J, GENERALCONSTANT) computes the kernel matrix for
%	the SDLFM kernel function in the block specified at indeces i,j. It
%	assumes the computation for functions that describe positions
%	(system 1) and velocity (system 2)
%	 Returns:
%	  K - the kernel matrix portion of block i,j
%	 Arguments:
%	  LFMKERN1 - structure containing parameters for the system 1
%	  LFMKERN2 - structure containing parameters for the system 2
%	  T1 - times at which the system 1 is evaluated
%	  T2 - times at which the system 2 is evaluated
%	  KYY - covariance for the initial conditions between position 1 and
%	   position 2 at block i,j
%	  KYV - covariance for the initial conditions between position 1 and
%	   velocity 2 at block i,j
%	  KVY - covariance for the initial conditions between velocity 1 and
%	   position 2 at block i,j
%	  KVV - covariance for the initial conditions between velocity 1 and
%	   velocity 2 at block i,j
%	  I - interval to be evaluated for system 1
%	  J - interval to be evaluated for system 2
%	  GENERALCONSTANT - constants evaluated with
%	   sdlfmKernComputeConstant.m


%	Copyright (c) 2010. Mauricio A. Alvarez


if nargin<11
    j = i;
    generalConst = [];
end

c1 = sdlfmMeanCompute(lfmKern1(1), t1, 'Pos');
e1 = sdlfmMeanCompute(lfmKern1(1), t1, 'Vel');
g2 = sdlfmvMeanCompute(lfmKern2(1), t2, 'Pos');
h2 = sdlfmvMeanCompute(lfmKern2(1), t2, 'Vel');

K = kyy*c1*g2.' + kyv*c1*h2.' + kvy*e1*g2.' + kvv*e1*h2.';

if i==j
    for k=1:length(lfmKern1)
        K  = K + (lfmvXlfmKernCompute(lfmKern2(k), lfmKern1(k), t2, t1).');
    end
else    
    if i>j
        PosVel = zeros(1,length(t2));
        VelVel = zeros(1,length(t2));
        for k=1:length(lfmKern1)
            PosVel = PosVel + lfmvXlfmKernCompute(lfmKern2(k), lfmKern1(k), t2, lfmKern2(k).limit).';
            VelVel = VelVel + lfmvXlfmvKernCompute(lfmKern1(k), lfmKern2(k), lfmKern2(k).limit, t2);
        end
        if isempty(generalConst{i,j})
            K = K + c1*PosVel+  e1*VelVel;
        else
            K = K + (generalConst{i,j}(1,1)*c1 + generalConst{i,j}(2,1)*e1)*PosVel + ...
                (generalConst{i,j}(1,2)*c1 + generalConst{i,j}(2,2)*e1)*VelVel;           
        end 
    else
        PosVel = zeros(length(t1),1);
        VelVel = zeros(length(t1),1);
        for k=1:length(lfmKern1)
            PosVel = PosVel + lfmXlfmKernCompute(lfmKern1(k), lfmKern2(k), t1, lfmKern1(k).limit);
            VelVel = VelVel + lfmvXlfmKernCompute(lfmKern2(k), lfmKern1(k), lfmKern1(k).limit, t1).';
            
        end        
         if isempty(generalConst{i,j})
             K = K + PosVel*g2.' + VelVel*h2.';               
         else
             K = K + PosVel*(generalConst{i,j}(1,1)*g2.' + generalConst{i,j}(2,1)*h2.') + ...
                 VelVel*(generalConst{i,j}(1,2)*g2.' + generalConst{i,j}(2,2)*h2.');
         end
    end
end

