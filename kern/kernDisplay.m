function kernDisplay(kern, varargin)

% KERNDISPLAY Display the parameters of the kernel.
% FORMAT
% DESC displays the parameters of the kernel and the kernel type to
% the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO modelDisplay
%
% COPYRIGHT Neil D. Lawrence, 2006, 2005, 2004

% KERN

fhandle = str2func([kern.type 'KernDisplay']);
fhandle(kern, varargin{:});
