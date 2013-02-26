function orderedNoiseDisplay(noise, spacing)

% ORDEREDNOISEDISPLAY Display parameters of the ORDERED noise.
%
%	Description:
%
%	ORDEREDNOISEDISPLAY(NOISE) displays the parameters of the ordered
%	categorical noise and the noise type to the console.
%	 Arguments:
%	  NOISE - the noise to display.
%
%	ORDEREDNOISEDISPLAY(NOISE, SPACING)
%	 Arguments:
%	  NOISE - the noise to display.
%	  SPACING - how many spaces to indent the display of the noise by.
%	
%
%	See also
%	ORDEREDNOISEPARAMINIT, MODELDISPLAY, NOISEDISPLAY


%	Copyright (c) 2004, 2005 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
for i = 1:noise.numProcess
  fprintf(spacing);
  fprintf('Ordered bias on process %d: %2.4f\n', i, noise.bias(i))
end
for i = 1:noise.C-2
  fprintf('Ordered noise model width %d: %2.4f\n', i, noise.widths(i))
end
fprintf(spacing);
fprintf(spacing);
fprintf('Ordered Sigma2: %2.4f\n', noise.variance);