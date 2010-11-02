 function g = multiKernDiagGradient(kern, x, covDiag)

% MULTIKERNDIAGGRADIENT Compute the gradient of the MULTI kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% multiple output block kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% multiKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% multiKernExtractParam.
%
% SEEALSO : multiKernParamInit, kernDiagGradient, multiKernExtractParam, multiKernGradient
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if iscell(x)
    dim1 = zeros(1, kern.numBlocks);
    % Collate arguments.
    for i=1:kern.numBlocks
        dim1(i) = size(x{i}, 1);
    end
    g = zeros(1, size(kern.paramGroups, 1));
    startVal = 1;
    endVal = 0;
    for i = 1:kern.numBlocks
        endVal = endVal + kern.comp{i}.nParams;
        startOne = sum(dim1(1:(i-1)))+1;
        endOne = sum(dim1(1:i));        
        g(1, startVal:endVal) = multiKernDiagGradientBlock(kern, x{i}, covDiag(startOne:endOne), i);
        startVal = endVal + 1;
    end
else
    % Collate arguments.
    dim1 = size(x, 1);
    arg{1} = x;
    g = zeros(1, size(kern.paramGroups, 1));
    startVal = 1;
    endVal = 0;
    for i = 1:kern.numBlocks
        endVal = endVal + kern.comp{i}.nParams;
        startOne = (i-1)*dim1 + 1;
        endOne = i*dim1;        
        g(1, startVal:endVal) = multiKernDiagGradientBlock(kern, arg{1}, covDiag(startOne:endOne), i);        
        startVal = endVal + 1;
    end

end

g = g*kern.paramGroups;


