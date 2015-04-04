function [kern, noise] = gpDeconstruct(model, fileName)

% GPDECONSTRUCT break GP in pieces for saving.

% GP

kern = rmfield(model.kern, 'Kstore');
kern = rmfield(kern, 'diagK');
noise = model.noise;