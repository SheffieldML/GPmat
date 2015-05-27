function [params, names] = sdrbfKernExtractParam(kern)

% SDRBFKERNEXTRACTPARAM Extract parameters from the SDRBF kernel structure.
% FORMAT
% DESC Extract parameters from the switching dynamical RBF kernel structure
% into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% FORMAT
% DESC Extract parameters and their names from the switching dynamical
% RBF kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO sdrbfKernParamInit, sdrbfKernExpandParam, kernExtractParam, 
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

params = [kern.inverseWidth(:)' kern.switchingTimes];
if nargout > 1
    namesInvWidth = cell(kern.nIntervals*kern.nlfPerInt,1);
    if kern.nlfPerInt == 1
        for i=1:kern.nIntervals
            namesInvWidth{i} = ['inverse width interval ' num2str(i) '.'];
        end
    else
        cont = 0;
        for i=1:kern.nlfPerInt
            for j=1:kern.nIntervals
                cont = cont + 1;
                namesInvWidth{cont} = ['inverse width ' num2str(i) '.' ' interval ' num2str(j) '.'];
            end
        end
    end
    namesStimes = cell(kern.nIntervals, 1);
    for i=1:kern.nIntervals
        namesStimes{i} = ['switching point interval ' num2str(i) '.'];
    end
    names = {namesInvWidth{:}, namesStimes{:}};
end
