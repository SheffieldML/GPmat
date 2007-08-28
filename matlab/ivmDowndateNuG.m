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

if isfield(model.kern, 'whiteVariance')
  subTerm = model.varSigma(index, :)-model.kern.whiteVariance;
else
  subTerm = model.varSigma(index, :);
end
subTerm = model.varSigma(index, :);
model.nu(index, :) = 1./(1./model.beta(index, :) - subTerm);
model.g(index, :) = model.nu(index, :).*(model.mu(index, :) - model.m(index, :));

