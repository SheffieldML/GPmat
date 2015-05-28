function sdlfmvKernDisplay(kern, spacing)

% SDLFMVKERNDISPLAY Display parameters of the SDLFMV kernel.
% FORMAT
% DESC displays the parameters of the switching dynamical latent force
% model kernel and the kernel type to the console. Displays the parameters
% for the velocity which are equal to the parameters for the position.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : sdlfmKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN

sdlfmKernDisplay(kern, spacing);
