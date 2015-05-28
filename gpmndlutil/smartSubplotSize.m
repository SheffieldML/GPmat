function [numRows, numColumns] = smartSubplotSize(N, row2colRatio)

% SMARTSUBPLOTSIZE 
% DESC Computes the smallest possible numbers numRows, numColumns to be used when calling
% subplot(numRows, numColumns, n) for the total number of plots being N.
% A preferance for the row/col ratio can be given optionally.
% COPYRIGHT: Andreas C. Damianou, 2013
% SHEFFIELDML

if nargin < 2
    row2colRatio = 3;
end

numRows = ceil(N.^(1/row2colRatio));
numColumns = ceil(N / numRows);
