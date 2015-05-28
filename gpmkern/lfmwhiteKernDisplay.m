function lfmwhiteKernDisplay(kern, spacing)

% LFMWHITEKERNDISPLAY Display parameters of the LFM-WHITE kernel.
% FORMAT
% DESC displays the parameters of the LFM-White (Latent Force Model - White)
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : lfmwhiteKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : David Luengo, 2009

% KERN


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
if kern.isStationary
    fprintf('Stationary version of the kernel\n');
else
    fprintf('Non-stationary version of the kernel\n');
end
fprintf(spacing);
fprintf('LFM-White mass: %2.4f\n', kern.mass);
fprintf(spacing);
fprintf('LFM-White spring: %2.4f\n', kern.spring);
fprintf(spacing);
fprintf('LFM-White damper: %2.4f\n', kern.damper);
fprintf(spacing);
fprintf('LFM-White Latent variance: %2.4f\n', kern.variance);
fprintf(spacing);
fprintf('LFM-White sensitivity: %2.4f\n', kern.sensitivity);
%fprintf(spacing);
%fprintf('LFM delay: %2.4f\n', kern.delay)
fprintf(spacing);
fprintf('System Characteristics:\n')
fprintf(spacing);
fprintf('LFM-White alpha: %2.4f\n', kern.alpha);
fprintf(spacing);
fprintf('LFM-White omega: %2.4f\n', kern.omega);
fprintf(spacing);
fprintf('LFM-White gamma: %2.4f\n', kern.gamma);
fprintf(spacing);
fprintf('LFM-White Damping Ratio: %2.4f\n', kern.zeta);
fprintf(spacing);
fprintf('LFM-White Undamped Natural Frequency: %2.4f\n', kern.omega_0);
