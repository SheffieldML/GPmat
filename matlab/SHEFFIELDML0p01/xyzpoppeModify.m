function handle = xyzpoppeModify(handle,pos)

% XYZPOPPEMODIFY
%
%	Description:
%	


%	Copyright (c) 2008 Carl Henrik Ek and Neil Lawrence



joint = xyzpoppe2joint(pos);

if(iscell(handle))
  for(i = 1:1:length(handle))
    xyzpoppeDraw(joint,handle{i});
  end
else
  xyzpoppeDraw(joint,handle);
end

return