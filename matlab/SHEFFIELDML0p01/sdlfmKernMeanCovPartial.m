function covGradLocal = sdlfmKernMeanCovPartial(lfmKern1, lfmKern2, t1, ...
    t2, covGrad, typeParam)

% SDLFMKERNMEANCOVPARTIAL Helper function for derivatives in SDLFM kernel
%
%	Description:
%	covGradLocal = sdlfmKernMeanCovPartial(lfmKern1, lfmKern2, t1, ...
%    t2, covGrad, typeParam)
%

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