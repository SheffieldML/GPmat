function [g1, g2, g3, g4] = sdlfmXsdlfmKernGradientIC(lfmKern1, lfmKern2, ...
    t1, t2, gkyy1, gkyy2, gkyy3, gkyy4, gkyv1, gkyv2, gkyv3, gkyv4, gkvy1, ...
    gkvy2, gkvy3, gkvy4, gkvv1, gkvv2, gkvv3, gkvv4, covGrad, typeParam)

% SDLFMXSDLFMKERNGRADIENTIC Computes partial derivative for init. const.
%
% COPYRIGHT : Mauricio A. Alvarez

% KERN

if nargin < 18
    typeParam = {'sdlfm', 'sdlfm'};    
end

g1 = cell(1,3);
g2 = cell(1,3);

covGradLocal = sdlfmKernMeanCovPartial(lfmKern1(1), lfmKern2(1), t1, t2, covGrad, ...
    typeParam);

% Parameters

for i=1:3,
    g1{i} = gkyy1{i}*covGradLocal(1,1) + gkyv1{i}*covGradLocal(1,2) ...
        + gkvy1{i}*covGradLocal(2,1) + gkvv1{i}*covGradLocal(2,2);
    g2{i} = gkyy2{i}*covGradLocal(1,1) + gkyv2{i}*covGradLocal(1,2) ...
        + gkvy2{i}*covGradLocal(2,1) + gkvv2{i}*covGradLocal(2,2);
end

% Switching points

g3 = gkyy3*covGradLocal(1,1) + gkyv3*covGradLocal(1,2) ...
    + gkvy3*covGradLocal(2,1) + gkvv3*covGradLocal(2,2);

% First covariance

g4 = gkyy4*covGradLocal(1,1) + gkyv4*covGradLocal(1,2) ...
   + gkvy4*covGradLocal(2,1) + gkvv4*covGradLocal(2,2);
