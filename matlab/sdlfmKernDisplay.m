function sdlfmKernDisplay(kern, spacing)

% SDLFMKERNDISPLAY Display parameters of the SDLFM kernel.
% FORMAT
% DESC displays the parameters of the switching dynamical latent force
% model kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : sdlfmKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010.

% KERN

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('SDLFM mass: %2.4f\n', kern.mass)
fprintf(spacing);
fprintf('SDLFM spring: %2.4f\n', kern.spring)
fprintf(spacing);
fprintf('SDLFM damper: %2.4f\n', kern.damper)
for i=1:kern.nlfPerInt
    for j=1:kern.nIntervals
        fprintf(spacing);
        fprintf('SDLFM inverse width %d interval %d: %2.4f (length scale %2.4f)\n', ...
           i,j, kern.inverseWidth(i,j), 1/sqrt(kern.inverseWidth(i,j)));        
    end
end
for j=1:kern.nIntervals
    fprintf(spacing);
    fprintf('SDLFM switching point interval %d: %2.4f \n', ...
        j, kern.switchingTimes(j));
end
for i=1:kern.nlfPerInt
    for j=1:kern.nIntervals
        fprintf(spacing);
        fprintf('SDLFM sensitivity %d interval %d: %2.4f \n', ...
           i,j, kern.sensitivity(i,j));        
    end
end
fprintf(spacing);
fprintf('System Characteristics:\n')
fprintf(spacing);
fprintf('SDLFM omega: %2.4f\n', kern.omega)
fprintf(spacing);
fprintf('SDLFM alpha: %2.4f\n', kern.alpha)
fprintf(spacing);
fprintf('SDLFM gamma: %2.4f\n', kern.gamma)
fprintf(spacing);
fprintf('SDLFM Damping Ratio: %2.4f\n', kern.zeta)
fprintf(spacing);
fprintf('SDLFM Undamped Natural Frequency: %2.4f\n', kern.omega_0)
