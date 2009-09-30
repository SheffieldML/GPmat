function [kern, noise, gpInfo] = gpDeconstruct(model)

% GPDECONSTRUCT break GP in pieces for saving.
% FORMAT
% DESC takes an GP model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN kern : the kernel component of the GP model.
% RETURN noise : the noise component of the GP model.
% RETURN gpInfo : a structure containing the other information
% from the GP: what the sparse approximation is, what the inducing
% variables are.
%
% SEEALSO : gpReconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009

% GP

kern = model.kern;
if isfield(model, 'noise')
  noise = model.noise;
else
  noise = [];
end
gpInfo = model;
removeFields = {'S', 'X', 'y', 'm', 'diagK', 'K_uu', 'invK_uu', 'logDetK_uu', ...
                'alpha', 'K_uf', 'sqrtK_uu', 'K', 'sqrtK_uu', 'innerProducts', ...
                'A', 'Ainv', 'logDetA', 'L', 'diagD', 'detDiff', 'Dinv', ...
                'logDetD', 'V', 'Am', 'Lm', 'invLmV', 'scaledM', 'bet', ...
                'D', 'noise', 'kern', 'expectations'};

for i = 1:length(removeFields)
  if isfield(gpInfo, removeFields{i})
    gpInfo = rmfield(gpInfo, removeFields{i});
  end
end

