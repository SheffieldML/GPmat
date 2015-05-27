function covGradLocal = sdlfmKernMeanCovPartial(lfmKern1, lfmKern2, t1, ...
    t2, covGrad, typeParam)

% SDLFMKERNMEANCOVPARTIAL Helper function for derivatives in SDLFM kernel
% FORMAT
% DESC computes partial derivatives for the total gradinets of the
% parameters in the SDLFM kernel
% ARG lfmKern1 : structure containing parameters system 1
% ARG lfmKern2 : structure containing parameters system 2
% ARG t1 : time points for the evaluation of system 1
% ARG t2 : time points for the evaluation of system 2
% ARG covGrad : portion of the partial derivative with respect the
% objective function that it's being optimised.
% ARG typeParam : associated to the type of kernel being computed.
% RETURN covGradLocal : the partial derivatives of the parameters

% KERN

if nargin < 6
    typeParam = {'sdlfm', 'sdlfm'};
end

fhandle1 = str2func([typeParam{1} 'MeanCompute']);
fhandle2 = str2func([typeParam{2} 'MeanCompute']);

c1 = fhandle1(lfmKern1, t1, 'Pos');
e1 = fhandle1(lfmKern1, t1, 'Vel');
c2 = fhandle2(lfmKern2, t2, 'Pos');
e2 = fhandle2(lfmKern2, t2, 'Vel');

covGradLocal(1,1) = sum(sum((c1*c2').*covGrad));
covGradLocal(1,2) = sum(sum((c1*e2').*covGrad));
covGradLocal(2,1) = sum(sum((e1*c2').*covGrad));
covGradLocal(2,2) = sum(sum((e1*e2').*covGrad));
