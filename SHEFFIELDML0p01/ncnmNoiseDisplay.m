function ncnmNoiseDisplay(noise, spacing)

% NCNMNOISEDISPLAY Display  parameters from null category noise model.
%
%	Description:
%
%	NCNMNOISEDISPLAY(NOISE) displays the parameters of the null category
%	noise model.
%	 Arguments:
%	  NOISE - the noise model to display.
%
%	NCNMNOISEDISPLAY(NOISE, SPACING)
%	 Arguments:
%	  NOISE - the noise model to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	% SEEALSO NCNMNOISEPARAMINIT, MODELDISPLAY, NOISEDISPLAY


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
for i = 1:noise.numProcess
  fprintf(spacing);
  fprintf('Null Cat bias on process %d: %2.4f\n', i, noise.bias(i))
end
fprintf(spacing);
fprintf('Null cat width: %2.4f\n', noise.width)
fprintf(spacing);
fprintf('Null cat Gamma -: %2.4f\n', noise.gamman)
fprintf(spacing);
fprintf('Null cat Gamma +: %2.4f\n', noise.gammap)
fprintf(spacing);
fprintf('Null cat Sigma2: %2.4f\n', noise.sigma2)
