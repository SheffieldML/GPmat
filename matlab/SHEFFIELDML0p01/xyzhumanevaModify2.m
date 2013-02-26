function handle = xyzhumanevaModify2(handle,pos)

% XYZHUMANEVAMODIFY2
%
%	Description:
%	


%	Copyright (c) 2008 Carl Henrik Ek and Neil Lawrence



joint = xyzhumaneva2joint(pos);
xyzhumanevaDraw(joint,handle{1});
xyzhumanevaDraw(joint,handle{2});

return