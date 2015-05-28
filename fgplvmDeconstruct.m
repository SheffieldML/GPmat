function [kern, noise, fgplvmInfo, X] = fgplvmDeconstruct(model)

% FGPLVMDECONSTRUCT break FGPLVM in pieces for saving.
% FORMAT
% DESC takes an FGPLVM model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN kern : the kernel component of the FGPLVM model.
% RETURN noise : the noise component of the FGPLVM model.
% RETURN fgplvmInfo : a structure containing the other information
% from the FGPLVM: what the active set is, what the inactive set is
% and what the site parameters are.
%
% SEEALSO : fgplvmReconstruct, gpDeconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2009

% FGPLVM

  [kern, noise, fgplvmInfo] = gpDeconstruct(model);
  if isfield(fgplvmInfo, 'back') && ~isempty(fgplvmInfo.back)
    switch fgplvmInfo.back.type
     case 'kbr'
      backToRemove = {'K', 'X'};
      
     otherwise
      backToRemove = {};
    end
    
    for i = 1:length(backToRemove)
      if isfield(fgplvmInfo.back, backToRemove{i})
        fgplvmInfo.back = rmfield(fgplvmInfo.back, backToRemove{i});
      end
    end
  end
  if isfield(fgplvmInfo, 'dynamics') && ~isempty(fgplvmInfo.dynamics)
    switch fgplvmInfo.dynamics.type 
     case 'gp'
      dynamicsToRemove = {};
     case 'gpDynamics'
      dynamicsToRemove = {'X', 'y', 'm', 'K_uu', 'K_uf', 'invK_uu', 'sqrtK_uu', 'logDetK_uu', ...
                          'diagK', 'L', 'diagD', 'Dinv', 'A', 'Ainv', 'logDetA', 'detDiff', ...
                          'innerProducts', 'V', 'Am', 'Lm', 'invLmV', ...
                          'scaledM', 'bet', 'K', 'D', 'logDetD'};
      
     otherwise
      dynamicsToRemove ={};
    end
    for i = 1:length(dynamicsToRemove)
      if isfield(fgplvmInfo.dynamics, dynamicsToRemove{i})
        fgplvmInfo.dynamics = rmfield(fgplvmInfo.dynamics, dynamicsToRemove{i});
      end
    end
  end
  X = model.X;
  fgplvmInfo.type = 'fgplvm';
end
