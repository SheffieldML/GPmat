function [vers, depend] = ivmVers

% IVMVERS Brings dependent toolboxes into the path.

vers = 0.32;
if nargout > 2
  depend(1).name = 'kern';
  depend(1).vers = 0.14;
  depend(1).required = 0;
  depend(2).name = 'noise';
  depend(2).vers = 0.13;
  depend(2).required = 0;
  depend(3).name = 'ndlutil';
  depend(3).vers = 0.13;
  depend(3).required = 0;
  depend(4).name = 'optimi';
  depend(4).vers = 0.12;
  depend(4).required = 0;
  depend(5).name = 'prior';
  depend(5).vers = 0.12;
  depend(5).required = 0;
  depend(6).name = 'rochol';
  depend(6).vers = 0.12;
  depend(6).required = 0;
end
