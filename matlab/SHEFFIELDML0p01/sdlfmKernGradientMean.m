function [grad1, grad2, gsp1, gsp2] = sdlfmKernGradientMean(lfmKern1, ...
    lfmKern2, t1, t2, kyy, kyv, kvy, kvv, covGrad, typeParam, typeSwitching)

% SDLFMKERNGRADIENTMEAN Gradients of the parameters of mean function cov.
%
%	Description:
%	[grad1, grad2, gsp1, gsp2] = sdlfmKernGradientMean(lfmKern1, ...
%    lfmKern2, t1, t2, kyy, kyv, kvy, kvv, covGrad, typeParam, typeSwitching)
%

if nargin<10
    typeParam{1} = 'sdlfm';
    typeParam{2} = 'sdlfm';
    typeSwitching{1} = 'sdlfmv';
    typeSwitching{2} = 'sdlfmv';
end

% Organizes the parameters

m = [lfmKern1.mass lfmKern2.mass];                  % Par. 1
D = [lfmKern1.spring lfmKern2.spring];              % Par. 2
C = [lfmKern1.damper lfmKern2.damper];              % Par. 3

alpha = C./(2*m);
omega = sqrt(D./m-alpha.^2);

% Derivatives of alpha and omega wrt parameters

gradAlphaParams1 = [-C(1)/(2*m(1)^2) 0 1/(2*m(1))];
gradAlphaParams2 = [-C(2)/(2*m(2)^2) 0 1/(2*m(2))];
gradOmegaParams1 = [(C(1)^2 - 2*m(1)*D(1))/(4*m(1)^3*omega(1)) ...
    1/(2*m(1)*omega(1)) -C(1)/(4*m(1)^2*omega(1))];
gradOmegaParams2 = [(C(2)^2 - 2*m(2)*D(2))/(4*m(2)^3*omega(2)) ...
    1/(2*m(2)*omega(2)) -C(2)/(4*m(2)^2*omega(2))];

% Computes the derivatives with respect to the parameters in the mean
% function.

fhandle1 = str2func([typeParam{1} 'MeanCompute']);
fhandle2 = str2func([typeParam{2} 'MeanCompute']);

c1 = fhandle1(lfmKern1, t1, 'Pos');
e1 = fhandle1(lfmKern1, t1, 'Vel');
c2 = fhandle2(lfmKern2, t2, 'Pos');
e2 = fhandle2(lfmKern2, t2, 'Vel');

fhandle3 = str2func([typeParam{1} 'MeanGradient']);
fhandle4 = str2func([typeParam{2} 'MeanGradient']);

[gc1Alpha, gc1Omega] = fhandle3(lfmKern1, t1, 'Pos');
[ge1Alpha, ge1Omega] = fhandle3(lfmKern1, t1, 'Vel');
[gc2Alpha, gc2Omega] = fhandle4(lfmKern2, t2, 'Pos');
[ge2Alpha, ge2Omega] = fhandle4(lfmKern2, t2, 'Vel');

% Derivative with respect to Alpha1 and Omega1

matGradAlpha = kyy*gc1Alpha*c2.' + kyv*gc1Alpha*e2.' + kvy*ge1Alpha*c2.' + ...
    kvv*ge1Alpha*e2.';

matGradOmega = kyy*gc1Omega*c2.' + kyv*gc1Omega*e2.' + kvy*ge1Omega*c2.' + ...
    kvv*ge1Omega*e2.';

grad1 = gradAlphaParams1*(sum(sum(matGradAlpha.*covGrad))) + ...
    gradOmegaParams1*(sum(sum(matGradOmega.*covGrad)));

% Derivatives with respect to Alpha2 and Omega2

matGradAlpha = kyy*c1*gc2Alpha.' + kyv*c1*ge2Alpha.' + kvy*e1*gc2Alpha.' + ...
    kvv*e1*ge2Alpha.';

matGradOmega = kyy*c1*gc2Omega.' + kyv*c1*ge2Omega.' + kvy*e1*gc2Omega.' + ...
    kvv*e1*ge2Omega.';

grad2 = gradAlphaParams2*(sum(sum(matGradAlpha.*covGrad))) + ...
    gradOmegaParams2*(sum(sum(matGradOmega.*covGrad)));

% Precomputations for the derivatives of the swicthing points

fhandle5 = str2func([typeSwitching{1} 'MeanCompute']);
fhandle6 = str2func([typeSwitching{2} 'MeanCompute']);

g1 = fhandle5(lfmKern1, t1, 'Pos');
h1 = fhandle5(lfmKern1, t1, 'Vel');
g2 = fhandle6(lfmKern2, t2, 'Pos');
h2 = fhandle6(lfmKern2, t2, 'Vel');

% Derivative of switching point 1

matGradSp1 = kyy*g1*c2.' + kyv*g1*e2.' + kvy*h1*c2.' + kvv*h1*e2.';

gsp1 = - sum(sum(matGradSp1.*covGrad));

% Derivative of switching point 2

matGradSp2 = kyy*c1*g2.' + kyv*c1*h2.' + kvy*e1*g2.' + kvv*e1*h2.';

gsp2 = - sum(sum(matGradSp2.*covGrad));

