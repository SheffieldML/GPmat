function ind = treeFindRoots(tree)

% TREEFINDROOTS Return indices of all root nodes in a tree structure.
% FORMAT
% DESC returns indices of all root nodes in an tree array.
% ARG tree : tree for which root nodes are being sought.
% RETURN ind : indices of root nodes.
%
% SEEALSO : treeFindParents, treeFindChildren, treeFindLeaves
%
% COPYRIGHT : Neil D. Lawrence, 2007

% NDLUTIL

ind = [];
for i = 1:length(tree)
  if isempty(tree(i).parent)
    ind = [ind i];
  end
end
