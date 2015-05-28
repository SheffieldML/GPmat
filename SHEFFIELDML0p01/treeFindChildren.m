function tree = treeFindChildren(tree)

% TREEFINDCHILDREN Given a tree that lists only parents, add children.
%
%	Description:
%
%	TREE = TREEFINDCHILDREN(TREE) takes a tree structure which lists the
%	children of each node and computes the parents for each node and
%	places them in.
%	 Returns:
%	  TREE - a tree that lists children and parents.
%	 Arguments:
%	  TREE - the tree that lists only children.
%	
%
%	See also
%	TREEFINDPARENTS


%	Copyright (c) 2005, 2006 Neil D. Lawrence


for i = 1:length(tree)
  for j = 1:length(tree(i).parent)
    if tree(i).parent(j)
      tree(tree(i).parent(j)).children ...
          = [tree(tree(i).parent(j)).children i];
    end
  end
end

