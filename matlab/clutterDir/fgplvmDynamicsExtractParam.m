function param = fgplvmDynamicsExtractParam(model)
 
% FGPLVMDYNAMICSEXTRACTPARAM Extract parameters from the dynamics portion.

fhandle = str2func([model.type 'ExtractParam']);

param = fhandle(model);

