function model = ivmDowndateNuG(model, index)

% IVMDOWNDATENUG Downdate nu and g parameters associated with noise model.

% IVM


model.nu(index, :) = 1./(1./model.beta(index, :) - model.varSigma(index, :));
model.g(index, :) = model.nu(index, :).*(model.mu(index, :) - model.m(index, :));

