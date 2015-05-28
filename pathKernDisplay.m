function pathKernDisplay(kern, spacing)

% PATHKERNDISPLAY Display parameters of the PATH kernel.
% FORMAT
% DESC displays the parameters of the path
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : pathKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek

% SHEFFIELDML


if nargin>1
  spacing=repmat(32,1,spacing);
else
  spacing=[];
end
spacing=char(spacing);
fprintf(spacing);
fprintf('Diagonal Step Cost: %2.4f\n',kern.cd);
fprintf(spacing);
fprintf('Hor/Vert Step Cost: %2.4f\n',kern.chv);
fprintf(spacing);
fprintf('Ground Kernel: %s\n',kern.gkern.type);
kernDisplay(kern.gkern,length(spacing)+4);
