function [kern, noise, ivmInfo] = ivmDeconstruct(model, fileName)

% IVMSAVE break IVM in pieces for saving.

% IVM

kern = rmfield(model.kern, 'Kstore');
kern = rmfield(kern, 'diagK');
noise = model.noise;
ivmInfo.I = model.I;
ivmInfo.J = model.J;
ivmInfo.m = model.m;
ivmInfo.beta = model.beta;

