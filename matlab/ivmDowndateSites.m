function model = ivmDowndateSites(model, index)

% IVMDOWNDATESITES Downdate site parameters.

% IVM

model.m(index, :) = 0;
model.beta(index, :) = 0;