function ngaussNoiseDisplay(noise)

% NGAUSSNOISEDISPLAY Display parameters of the NGAUSS noise.
%
%	Description:
%
%	NGAUSSNOISEDISPLAY(NOISE) displays the parameters of the noiseless
%	Gaussian noise and the noise type to the console.
%	 Arguments:
%	  NOISE - the noise to display.
%
%	NGAUSSNOISEDISPLAY(NOISE, SPACING)
%	 Arguments:
%	  NOISE - the noise to display.
%	  SPACING - how many spaces to indent the display of the noise by.
%	
%
%	See also
%	NGAUSSNOISEPARAMINIT, MODELDISPLAY, NOISEDISPLAY


%	Copyright (c) 2004, 2005 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
for i = 1:noise.numProcess
  fprintf(spacing)
  fprintf('Noiseless Gaussian bias on process %d: %2.4f\n', i, noise.bias(i))
end
