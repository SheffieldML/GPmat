function g = ouKernGradient(kern, t1, varargin)

% OUKERNGRADIENT Gradient of OU kernel's parameters (see ouKernCompute or
%
%	Description:
%	ouKernParamInit for a more detailed description of the OU kernel).
%
%	G = OUKERNGRADIENT(KERN, T1, PARTIAL) computes the gradient of
%	functions with respect to the Ornstein-Uhlenbeck kernel's
%	parameters. As well as the kernel structure and the input positions,
%	the user provides a matrix PARTIAL which gives the partial
%	derivatives of the function with respect to the relevant elements of
%	the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  T1 - the input locations for which the gradients are being
%	   computed.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in t1.
%
%	G = OUKERNGRADIENT(KERN, T1, T2, PARTIAL) computes the derivatives
%	as above, but input locations are now provided in two matrices
%	associated with rows and columns of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  T1 - the input locations associated with the rows of the kernel
%	   matrix.
%	  T2 - the input locations associated with the columns of the kernel
%	   matrix.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as t1 and the same number of columns as t2
%	   has rows.
%	
%
%	See also
%	% SEEALSO OUKERNPARAMINIT, KERNGRADIENT, OUKERNDIAGGRADIENT, KERNGRADX


%	Copyright (c) 2009 David Luengo



if nargin < 4
    t2 = t1;
else
    t2 = varargin{1};
end

if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
K1 = exp(-kern.decay*abs(T1-T2));
if (kern.isStationary == false)
    K2 = exp(-kern.decay*(T1+T2));
else
    K2 = 0;
end
g(1) = (0.5*kern.variance/kern.decay)*sum(sum(varargin{end}.* ...
    ((K2-K1)/kern.decay - abs(T1-T2).*K1 + (T1+T2).*K2)));
g(2) = (0.5/kern.decay)*sum(sum(varargin{end}.*(K1-K2)));
