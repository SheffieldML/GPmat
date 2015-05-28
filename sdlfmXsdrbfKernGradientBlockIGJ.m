function [g1, g2, g3] = sdlfmXsdrbfKernGradientBlockIGJ(lfmKern, rbfKern, ...
    t1, t2, i, j, generalConst, generalConstGrad, covGrad, typeMeanGrad)

% SDLFMXSDRBFKERNGRADIENTBLOCKIGJ 
% FORMAT
% DESC computes the gradients of the parameters for output system and
% latent force system when i is greater than j.
% ARG lfmKern : structure containing parameters for the output system
% ARG rbfKern : structure containing parameters for the latent force system 
% ARG t1 : times at which the output system is evaluated
% ARG t2 : times at which the LF system is evaluated
% ARG i : interval in which the output is evaluated
% ARG j : interval in which the latent force is being evaluated
% ARG i : interval to be evaluated for system 1
% ARG j : interval to be evaluated for system 2
% ARG generalConstant : constants evaluated with sdlfmKernComputeConstant.m
% ARG generalConstGrad : derivatives of the constants computed with
% sdlfmKernGradientConstant.m
% ARG covGrad : partial derivatives of the objective function wrt portion
% of the corresponding kernel matrix
% ARG typeMeanGrad : specify the mean functions used to compute this part 
% of the kernel
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN

if nargin < 10
    typeMeanGrad = {'lfm', 'lfmv'};
end

fhandleMeanParam = str2func(['sd' typeMeanGrad{1} 'MeanCompute']);
fhandleMeanGradParam = str2func(['sd' typeMeanGrad{1} 'MeanGradient']);
fhandleMeanGradSwitching = str2func(['sd' typeMeanGrad{2} 'MeanCompute']);

c1 = fhandleMeanParam(lfmKern, t1, 'Pos');
e1 = fhandleMeanParam(lfmKern, t1, 'Vel');

if isempty(generalConst{i,j})
    coeffPosRbf = c1;
    coeffVelRbf = e1;
else
    coeffPosRbf = generalConst{i,j}(1,1)*c1 + generalConst{i,j}(2,1)*e1;
    coeffVelRbf = generalConst{i,j}(1,2)*c1 + generalConst{i,j}(2,2)*e1;
end

PosRbf = lfmXrbfKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
[g1Kern, g2Kern] = lfmXrbfKernGradient(lfmKern, rbfKern, ...
        rbfKern.limit, t2, covGrad.', coeffPosRbf);
VelRbf = lfmvXrbfKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
[g1KLocal, g2KLocal] = lfmvXrbfKernGradient(lfmKern, rbfKern, ...
        rbfKern.limit, t2, covGrad.', coeffVelRbf);
g1Kern = g1Kern + g1KLocal;
g2Kern = g2Kern + g2KLocal;

[gcAlpha, gcOmega] = fhandleMeanGradParam(lfmKern, t1, 'Pos');
[geAlpha, geOmega] = fhandleMeanGradParam(lfmKern, t1, 'Vel');
[gradAlpha, gradOmega] = getLocalGradAlphaOmega(lfmKern);

g1Pos = fhandleMeanGradSwitching(lfmKern, t1, 'Pos');
h1Vel = fhandleMeanGradSwitching(lfmKern, t1, 'Vel');

% This means that are only two switching points involved, t1
% and t0 appear in k, like k_{ff}(t1-t0,t'-t0) and
% k_{mf}(t1-t0,t'-t0). The other points appear in c1(t - t1)
% and e1(t - t1)
gsp1Kern = 0;
gsp2Kern = 0;

temp = lfmvXrbfKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
covPartial = sum(sum((coeffPosRbf*temp).*covGrad));
gsp2Kern = gsp2Kern - covPartial;
gsp1Kern = gsp1Kern + covPartial;
temp = lfmXrbfvKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
covPartial = sum(sum((coeffPosRbf*temp).*covGrad));
gsp2Kern = gsp2Kern - covPartial;
temp = lfmaXrbfKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
covPartial = sum(sum((coeffVelRbf*temp).*covGrad));
gsp2Kern = gsp2Kern - covPartial; 
gsp1Kern = gsp1Kern + covPartial;
temp = lfmvXrbfvKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
gsp2Kern = gsp2Kern - sum(sum((coeffVelRbf*temp).*covGrad));

% Organize gradients accordingly
if isempty(generalConst{i,j})
    matGradAlpha = gcAlpha*PosRbf + geAlpha*VelRbf;
    matGradOmega = gcOmega*PosRbf + geOmega*VelRbf;
    gCoeff = gradAlpha*sum(sum(matGradAlpha.*covGrad)) + ...
        gradOmega*(sum(sum(matGradOmega.*covGrad)));
    matGradSp1 = g1Pos*PosRbf + h1Vel*VelRbf;
    gsp1 = - sum(sum(matGradSp1.*covGrad));
    g3(i) = sum(gsp1Kern) + gsp1; % switching point 1
    g3(j) = sum(gsp2Kern);        % switching point 2
else
    constGradAlpha = generalConstGrad{1}{i,j};
    constGradOmega = generalConstGrad{2}{i,j};
    constVal = generalConst{i,j};
    % Derivative wrt parameters
    matGradAlpha = (constVal(1,1)*gcAlpha + constGradAlpha(1,1)*c1 ...
        + constVal(2,1)*geAlpha + constGradAlpha(2,1)*e1)*PosRbf ...
        + (constVal(1,2)*gcAlpha + constGradAlpha(1,2)*c1 ...
        + constVal(2,2)*geAlpha + constGradAlpha(2,2)*e1)*VelRbf;
    matGradOmega = (constVal(1,1)*gcOmega + constGradOmega(1,1)*c1 ...
        + constVal(2,1)*geOmega + constGradOmega(2,1)*e1)*PosRbf ...
        + (constVal(1,2)*gcOmega + constGradOmega(1,2)*c1 ...
        + constVal(2,2)*geOmega + constGradOmega(2,2)*e1)*VelRbf;
    gCoeff = gradAlpha*sum(sum(matGradAlpha.*covGrad)) + ...
        gradOmega*(sum(sum(matGradOmega.*covGrad)));
    % Derivative wrt switching points
    % Firts, notice that gsp1Kern doesn't correspond to the switching point
    % of the interval i (or say 1). It's just the derivative of the inner
    % switching point that leads to the actual switching point for the
    % interval 1. So we rename it first.
    gspInt = sum(gsp1Kern);
    % Compute the derivative of the swicthing point i, not appearing in the
    % constant
    matGradSp1 = (constVal(1,1)*g1Pos + constVal(2,1)*h1Vel)*PosRbf + ...
        (constVal(1,2)*g1Pos + constVal(2,2)*h1Vel)*VelRbf;
    gsp1 = - sum(sum(matGradSp1.*covGrad));
    % Compute the derivatives for all the switching points appearing in the
    % constant. The order is descending, i.e., t_3_, t_2, t_1
    constGradSPoint = generalConstGrad{3}{i,j};
    numberSP = size(constGradSPoint,2);
    gspInBetween = zeros(1, numberSP);
    for k=1:numberSP
        temp = (constGradSPoint(1,k)*c1 + constGradSPoint(3,k)*e1)*PosRbf + ...
            (constGradSPoint(2,k)*c1 + constGradSPoint(4,k)*e1)*VelRbf;
        gspInBetween(k) = sum(sum(temp.*covGrad));
    end
    % Assign derivatives wrt all other switching points, with a correction
    % for the derivative of the innermost switching point in the constant
    gspInBetween(end) =  gspInBetween(end) + gspInt;
    g3(i) = gsp1 + gspInBetween(1);
    g3(j) = sum(gsp2Kern);
    g3(j+1:i-1) = fliplr(gspInBetween(2:end));
end
g1 = zeros(1,5);
g2 = zeros(1,2);
% Assign derivatives wrt first system
g1(1:3) = g1Kern(1:3) + gCoeff;          % mass, spring, damper
g1(4) = g1Kern(4);                            % inverse width
g1(5) = g1Kern(5);                            % sensitivity
% Assign derivatives wrt second system
g2(1) = g2Kern(1);                      % variance
g2(2) = g2Kern(2);                   % inverse widths
