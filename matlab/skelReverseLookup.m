function index = skelReverseLookup(skel, jointName)

% SKELREVERSELOOKUP Return the number associated with the joint name.

% MOCAP

for i=1:length(skel.tree)
  if strcmp(skel.tree(i).name, jointName)
    index = i;
    return
  end
end

error('Reverse look up of name failed.')
