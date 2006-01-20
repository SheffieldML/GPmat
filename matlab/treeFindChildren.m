function tree = treeFindChildren(tree)

% TREEFINDCHILDREN Given a tree that lists only parents, add children.

% MOCAP

for i = 1:length(tree)
  for j = 1:length(tree(i).parent)
    if tree(i).parent(j)
      tree(tree(i).parent(j)).children ...
          = [tree(tree(i).parent(j)).children i];
    end
  end
end

