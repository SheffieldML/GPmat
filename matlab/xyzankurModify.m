function handle = xyzankurModify(handle,pos)

% XYZANKURMODIFY
%
% COPYRIGHT : Carl Henrik Ek and Neil Lawrence, 2008

% VISUALISATION


joint = xyzankur2joint(pos);

if(iscell(handle))
  for(i = 1:1:length(handle))
    xyzankurDraw(joint,handle{i}); 
  end
else
  xyzankurDraw(joint,handle);
end

return