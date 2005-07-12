function kern = kernExpandParam(kern, params)

% KERNEXPANDPARAM Expand parameters to form a kernel structure.

% KERN

% Check if parameters are being optimised in a transformed space.
if ~isempty(kern.transforms)
  for i = 1:length(kern.transforms)
    index = kern.transforms(i).index;
    fhandle = str2func([kern.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'atox');
  end
end
fhandle = str2func([kern.type 'KernExpandParam']);
kern = fhandle(kern, params);

