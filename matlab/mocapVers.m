function [vers, depend] = mocapVers

% MOCAPPATH Brings dependent toolboxes into the path.

% MOCAP

vers = 0.001;
if nargout > 2
  depend(1).name = 'gplvm';
  depend(1).vers = 1.01;
  depend(1).required = 1;
end
