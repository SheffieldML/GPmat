function [x, t] = timeseriesdata(data, nin);

% TIMESERIESDATA make a time series data set with the given window length.

% NDLUTIL

ndata = length(data);
alldata = ones(ndata-nin+1, nin+1); % [y_{t-1], ..., y_{t-m]]

% Creates the delay vectors
for j = 1:ndata-nin
  alldata(j, :) = (data(j:j+nin)');
end

% Training set
x = alldata(:,1:nin);
t = alldata(:,nin+1);
