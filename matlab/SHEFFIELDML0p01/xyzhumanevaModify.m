function handle = xyzhumanevaModify(handle,pos)

% XYZHUMANEVAMODIFY
%
%	Description:
%	


%	Copyright (c) 2008 Carl Henrik Ek and Neil Lawrence



joint = xyzhumaneva2joint(pos);

if(iscell(handle))
  for(i = 1:1:length(handle))
    xyzhumanevaDraw(joint,handle{i});
  end
else
  xyzhumanevaDraw(joint,handle);
end

return