function tree = treeFindParents(tree)

% TREEFINDPARENTS Given a tree that lists only children, add parents.
% FORMAT 
% DESC takes a tree structure which lists the parents of each node
% and computes the children for each node and places them in.
% ARG tree : the tree that lists only parents.
% RETURN tree : a tree that lists parents and children.
%
% SEEALSO : treeFindChildren
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% NDLUTIL

for i = 1:length(tree)
  for j = 1:length(tree(i).children)
    if tree(i).children(j)
      tree(tree(i).children(j)).parent ...
          = [tree(tree(i).children(j)).parent i];
    end
  end
end

