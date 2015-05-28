function [params, names] = gaussian2PriorExtractParam(prior)

% GAUSSIAN2PRIOREXTRACTPARAM Extract params from Gaussian prior structure.
% COPYRIGHT: Neil D. Lawrence 2004
% MODIFICATIONS: Andreas C. Damianou, 2013
%
% SHEFFIELDML

params = prior.precision;
params = [params prior.mean'];

if nargout > 1
    names = {'Gaussian precision'};
    for i = 1:length(prior.mean)
        names{i+1} = ['Gaussian mean ' num2str(i)];
    end
end
