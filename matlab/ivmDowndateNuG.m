function model = ivmDowndateNuG(model, index)

% IVMDOWNDATENUG Downdate nu and g parameters associated with noise model.
% FORMAT
% DESC removes a given data point from the nu and g
% representations.
% ARG model : the model from which the point is to be removed.
% ARG index : the index of the point to be removed.
% RETURN model : the model with the point removed.
% 
% SEEALSO : ivmRemovePoint, ivmEpUpdatePoint, ivmDowndateM
%
% COPYRIGHT : Neil D. Lawrence, 2005

% IVM


model.nu(index, :) = 1./(1./model.beta(index, :) - model.varSigma(index, :));
model.g(index, :) = model.nu(index, :).*(model.mu(index, :) - model.m(index, :));

