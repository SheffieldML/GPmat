function g = kernGradient(kern, x, covGrad)

% KERNGRADIENT Compute the gradient of the kernel's parameters.

% IVM

g = feval([kern.type 'KernGradient'], kern, x, covGrad);

% Check if parameters are being optimised in a transformed space.
if isfield(kern, 'transforms')
  params = feval([kern.type 'KernExtractParam'], kern);
  for i = 1:length(kern.transforms)
    index = kern.transforms(i).index;
    g(index) = g(index).* ...
        feval([kern.transforms(i).type 'Transform'], ...
              params(index), 'gradfact');
  end
end

  