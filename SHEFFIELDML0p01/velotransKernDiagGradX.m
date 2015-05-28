function gX = velotransKernDiagGradX(kern, X)

% VELOTRANSKERNDIAGGRADX Gradient of VELOTRANS kernel's diagonal with respect to X.
%
%	Description:
%
%	GX = VELOTRANSKERNDIAGGRADX(KERN, X) computes the gradient of the
%	diagonal of the velocity translate kernel matrix with respect to the
%	elements of the design matrix given in X.
%	 Returns:
%	  GX - the gradients of the diagonal with respect to each element of
%	   X. The returned matrix has the same dimensions as X.
%	 Arguments:
%	  KERN - the kernel structure for which gradients are being
%	   computed.
%	  X - the input data in the form of a design matrix.
%	
%
%	See also
%	VELOTRANSKERNPARAMINIT, KERNDIAGGRADX, VELOTRANSKERNGRADX, TRANSLATEKERNDIAGGRADX


%	Copyright (c) 2011 Neil D. Lawrence


  Xpass = X(:, 1:end-1);
  t = X(:, end);
  Xpass = Xpass - t*kern.velocity;
  gX = cmpndKernDiagGradX(kern, Xpass);
  gX = [gX  -gX*kern.velocity']; % That's a wild guess at the gradient, but
                               % Geoff Hinton told me he guessed the t-SNE
                               % gradients so that's O.K. :-)
end