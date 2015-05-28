function out = cell2num(inCell)

% CELL2NUM Converts a cell array of strings to numbers.
% FORMAT
% DESC converts a cell array of strings to numbers.
% ARG inCell : the input cell array.
% RETURN out : the output array.
%
% SEEALSO : str2num, cell2mat
%
% COPYRIGHT : Neil D. Lawrence, 2010

% NDLUTIL

  if nargin < 3
    before = true;
    if nargin < 2
      padding = ' ';
    end
  end

  lengths = cellfun('length', inCell);
  maxLength = max(lengths);
  maxPad = maxLength - min(lengths);
  out = zeros(length(inCell), 1);
  for i = 0:maxPad
    ind = find(lengths==maxLength-i);
    if length(ind)==0
      continue
    end
    out(ind) = str2num(cat(1, inCell{ind}));
  end
end
