function [params, names] = sdlfmKernExtractParam(kern)

% SDLFMKERNEXTRACTPARAM Extract parameters from the SDLFM kernel structure.
% FORMAT
% DESC Extract parameters from the switching dynamical LFM kernel structure
% into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% FORMAT
% DESC Extract parameters and their names from the switching dynamical
% LFM kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO sdlfmKernParamInit, sdlfmKernExpandParam, kernExtractParam, 
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

params = [kern.mass, kern.spring, kern.damper, ...
    kern.inverseWidth(:)', kern.switchingTimes, kern.sensitivity(:)'];
if nargout > 1
    names = {'mass', 'spring', 'damper'};
    namesInvWidth = cell(kern.nIntervals*kern.nlfPerInt,1);
    namesSensitivities = cell(kern.nIntervals*kern.nlfPerInt,1);
    if kern.nlfPerInt == 1
        for i=1:kern.nIntervals
            namesInvWidth{i} = ['inverse width interval ' num2str(i) '.'];
            namesSensitivities{i} = ['sensitivity interval ' num2str(i) '.'];
        end
    else
        cont = 0;
        for i=1:kern.nlfPerInt
            for j=1:kern.nIntervals
                cont = cont + 1;
                namesInvWidth{cont} = ['inverse width ' num2str(i) '.' ' interval ' num2str(j) '.'];
                namesSensitivities{cont} = ['sensitivity ' num2str(i) '.' ' interval ' num2str(j) '.'];
            end
        end
    end
    namesStimes = cell(kern.nIntervals, 1);
    for i=1:kern.nIntervals
        namesStimes{i} = ['switching point interval ' num2str(i) '.'];
    end
    names = {names{:}, namesInvWidth{:}, namesStimes{:}, ...
        namesSensitivities{:}};
end
