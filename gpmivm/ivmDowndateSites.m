function model = ivmDowndateSites(model, index)

% IVMDOWNDATESITES Downdate site parameters.
% FORMAT
% DESC remove a given data point from the site parameters.
% ARG model : the IVM structure from which the point is to be
% removed.
% ARG index : the index of the point to be removed.
% RETURN model : the model with the given point removed.
%
% COPYRIGHT : Neil D. Lawrence, 2005
%
% SEEALSO : ivmEpUpdateM, ivmRemovePoint

% IVM

model.m(index, :) = 0;
model.beta(index, :) = 0;
