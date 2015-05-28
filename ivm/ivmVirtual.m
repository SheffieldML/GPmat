function [X, y] = ivmVirtual(origX, origy, invariance)

% IVMVIRTUAL Create virtual data points with the specified invariance.
% One advantage of a sparse compression scheme (like the IVM and
% SVM, but unlike the RVM) is you can impose invariances on your
% data by only placing invariances on the 'active set'. This is a
% simple helper function to create a new data set (typically from
% an active set) by imposing a simple invariance.
%
% FORMAT
% DESC creates virtual data points for use with a 'virtual IVM'. 
% ARG X : input data set on which we wish to make the transformation
% IVM

switch invariance
 case 'translate'
  baseVal = min(min(origX));
  X = zeros(size(origX, 1)*5, size(origX, 2));
  y = zeros(size(origy, 1)*5, size(origy, 2));
  for j = 1:size(origX, 1)
    xPoint =  origX(j, :);
    X(5*j-4, :) = xPoint;
    Xmat = reshape(xPoint, 16, 16);
    XmatU = [Xmat(2:end, :); repmat(baseVal, 1, 16)];
    X(5*j-3, :) = XmatU(:)';
    XmatD = [repmat(baseVal, 1, 16); Xmat(1:end-1, :)];
    X(5*j-2, :) = XmatD(:)';
    XmatL = [Xmat(:, 2:end) repmat(baseVal, 16, 1)];
    X(5*j-1, :) = XmatL(:)';
    XmatR = [repmat(baseVal, 16, 1) Xmat(:, 1:end-1)];
    X(5*j, :) = XmatR(:)';
    for i = 0:4
      y(5*j-i, :) = origy(j, :);
    end
  end
 otherwise
  error('That invariance is not yet implemented')
end
