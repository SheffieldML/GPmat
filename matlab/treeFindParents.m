function tree = treeFindParents(tree)

% TREEFINDPARENTS Given a tree that lists only children, add parents.

% MOCAP

for i = 1:length(tree)
  for j = 1:length(tree(i).children)
    if tree(i).children(j)
      tree(tree(i).children(j)).parent ...
          = [tree(tree(i).children(j)).parent i];
    end
  end
end

