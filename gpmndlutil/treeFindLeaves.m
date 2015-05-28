function ind = treeFindLeaves(tree)

% TREEFINDLEAVES Return indices of all leaf nodes in a tree structure.
% FORMAT
% DESC returns indices of all leaf nodes in an tree array.
% ARG tree : tree for which leaf nodes are being sought.
% RETURN ind : indices of leaf nodes.
%
% SEEALSO : treeFindParents, treeFindChildren
%
% COPYRIGHT : Neil D. Lawrence, 2007

% NDLUTIL

ind = [];
for i = 1:length(tree)
  if isempty(tree(i).children)
    ind = [ind i];
  end
end
