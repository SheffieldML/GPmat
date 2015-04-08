function [vers, depend] = kernVers

% KERNPATH Brings dependent toolboxes into the path.

vers = 0.131;
if nargout > 2
  depend(1).name = 'ndlutil';
  depend(1).vers = 0.13;
  depend(1).required = 0;
end
