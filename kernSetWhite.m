function kern = kernSetWhite(kern, value)

% KERNSETWHITE Helper function to set the white noise in a kernel if it exists.

% KERN

whitePresent = 0;
switch kern.type
 case 'cmpnd'
  if isfield(kern, 'whiteVariance')
    kern.whiteVariance = value;
    for i = 1:length(kern.comp)
      if strcmp(kern.comp{i}.type, 'white') 
        if whitePresent
          error(['There are multiple sources of white noise in kernel, ' ...
                 'cannot set white noise uniquely.']) 
        end
        kern.comp{i}.variance = value;
        whitePresent = 1;
      end
      if isfield(kern.comp{i}, 'whiteVariance')
        if whitePresent
          error(['There are multiple sources of white noise in kernel, ' ...
                 'cannot set white noise uniquely.'])
        end
        kern.comp{i}.whiteVariance = value;
      end
    end
  end
  
 case 'white'
  kern.variance = value;
  whitePresent = 1;
 
 otherwise
  if isfield(kern, 'whiteVariance')
    kern.whiteVariance = value;
    whitePresent = 1;
  end
end

if ~whitePresent
  warning(['Attempted to set white noise in kernel when no white noise ' ...
           'term is present'])
end
