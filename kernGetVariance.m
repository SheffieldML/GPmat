function var = kernGetVariance(kern)

% KERNGETVARIANCE Get the signal associated with a the kernel.

% KERN

fhandleStr = ([kern.type 'KernGetVariance']);
if exist(fhandleStr)==2
  fhandle = str2func(fhandleStr);
  var = fhandle(kern);
  return
else
  if(strcmp(kern.type, 'white'))
    var = 0;
    return
  else
    if(isfield(kern, 'variance'))
      var = kern.variance;
    else
      error(['No field variance or file ' fhandleStr])
    end
    return
  end
end
