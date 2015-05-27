function Kc = kernelCenter(K)

% KERNELCENTER Attempts to Center Kernel Matrix
% FORMAT
% DESC returns a centered kernel matrix
% ARG kernel matrix
% RETURN Kc : The centered kernel
%
% SEEALSO : 
%
% COPYRIGHT : Carl Henrik Ek, 2008

% KERN

if(nargin<1)
  error('To Few Arguments');
end
if(size(K,1)~=size(K,2))
  error('Kernel Not Square');
end

Kc = (eye(size(K)) - 1/size(K,1)*ones(size(K)))*K*(eye(size(K))-1/size(K,1)*ones(size(K)));

return;
