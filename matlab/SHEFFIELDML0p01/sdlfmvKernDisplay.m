function sdlfmvKernDisplay(kern, spacing)

% SDLFMVKERNDISPLAY Display parameters of the SDLFMV kernel.
%
%	Description:
%
%	SDLFMVKERNDISPLAY(KERN) displays the parameters of the switching
%	dynamical latent force model kernel and the kernel type to the
%	console. Displays the parameters for the velocity which are equal to
%	the parameters for the position.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	SDLFMVKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	SDLFMKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2010. Mauricio A. Alvarez


sdlfmKernDisplay(kern, spacing);