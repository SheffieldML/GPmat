function [vers, depend] = rocholVers

% ROCHOLPATH Brings dependent toolboxes into the path.

% ROCHOL

vers = 0.13;
if nargout > 2
  depend(1).name = 'ndlutil';
  depend(1).vers = 0.13;
  depend(1).required = 0;
end
