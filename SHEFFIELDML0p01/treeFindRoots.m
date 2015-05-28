function ind = treeFindRoots(tree)

% TREEFINDROOTS Return indices of all root nodes in a tree structure.
%
%	Description:
%
%	IND = TREEFINDROOTS(TREE) returns indices of all root nodes in an
%	tree array.
%	 Returns:
%	  IND - indices of root nodes.
%	 Arguments:
%	  TREE - tree for which root nodes are being sought.
%	
%
%	See also
%	TREEFINDPARENTS, TREEFINDCHILDREN, TREEFINDLEAVES


%	Copyright (c) 2007 Neil D. Lawrence


ind = [];
for i = 1:length(tree)
  if isempty(tree(i).parent)
    ind = [ind i];
  end
end