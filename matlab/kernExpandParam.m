function kern = kernExpandParam(kern, params)

% KERNEXPANDPARAM Expand parameters to form a kernel structure.

% IVM

% Check if parameters are being optimised in a transformed space.
if isfield(kern, 'transforms')
  for i = 1:length(kern.transforms)
    index = kern.transforms(i).index;
    params(index) = feval([kern.transforms(i).type 'Transform'], ...
              params(index), 'atox');
  end
end

kern = feval([kern.type 'KernExpandParam'], kern, params);

