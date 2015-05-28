function [kern, noise, ivmInfo] = ivmDeconstruct(model)

% IVMDECONSTRUCT break IVM in pieces for saving.
% FORMAT
% DESC takes an IVM model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN kern : the kernel component of the IVM model.
% RETURN noise : the noise component of the IVM model.
% RETURN ivmInfo : a structure containing the other information
% from the IVM: what the active set is, what the inactive set is
% and what the site parameters are.
%
% SEEALSO : ivmReconstruct
%
% COPYRIGHT : Neil D. Lawrence, 2005

% IVM

kern = rmfield(model.kern, 'Kstore');
kern = rmfield(kern, 'diagK');
noise = model.noise;
ivmInfo.I = model.I;
ivmInfo.J = model.J;
ivmInfo.m = model.m;
ivmInfo.beta = model.beta;
ivmInfo.type = 'ivm';
