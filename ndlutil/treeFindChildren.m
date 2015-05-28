function tree = treeFindChildren(tree)

% TREEFINDCHILDREN Given a tree that lists only parents, add children.
% FORMAT 
% DESC takes a tree structure which lists the children of each node
% and computes the parents for each node and places them in.
% ARG tree : the tree that lists only children.
% RETURN tree : a tree that lists children and parents.
%
% SEEALSO : treeFindParents
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% NDLUTIL

for i = 1:length(tree)
  for j = 1:length(tree(i).parent)
    if tree(i).parent(j)
      tree(tree(i).parent(j)).children ...
          = [tree(tree(i).parent(j)).children i];
    end
  end
end

