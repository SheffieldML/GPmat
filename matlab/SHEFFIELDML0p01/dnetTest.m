
% DNETTEST Test some settings for the density network.
%
%	Description:
%	

model = dnetCreate(2, 3, randn(2, 3), dnetOptions);

model.basisStored = false;

modelTest(model)  