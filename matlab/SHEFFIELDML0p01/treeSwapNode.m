function tree = treeSwapNode(tree, i, j);

% TREESWAPNODE Swap two nodes in the tree structure array.
%
%	Description:
%
%	TREE = TREESWAPNODE(TREE, I, J) swaps the location of two nodes in a
%	tree structure array.
%	 Returns:
%	  TREE - the tree structure with the two node locations swapped.
%	 Arguments:
%	  TREE - the tree for which two nodes are to be swapped.
%	  I - the index of the first node to be swapped.
%	  J - the index of the second node to be swapped.
%	
%
%	See also
%	TREEFINDPARENTS, TREEFINDCHILDREN


%	Copyright (c) 2005, 2006 Neil D. Lawrence


storeNodeI = tree(i);
storeNodeJ = tree(j);
tree(j) = storeNodeI;
tree(i) = storeNodeJ;
for k = 1:length(tree)
  tree(k).children(find(tree(k).children==i)) = -1;
  tree(k).children(find(tree(k).children==j)) = i;
  tree(k).children(find(tree(k).children==-1)) = j;
  tree(k).parent(find(tree(k).parent==i)) = -1;
  tree(k).parent(find(tree(k).parent==j)) = i;
  tree(k).parent(find(tree(k).parent==-1)) = j;
end
