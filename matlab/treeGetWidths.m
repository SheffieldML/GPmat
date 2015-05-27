function [widths, maxDepth, nodePositions] = treeGetWidths(tree)

% TREEGETWIDTHS give width of each level of tree.
% FORMAT
% DESC gives the width of a tree at each level of the hierarchy and
% the maximum depth of a tree as well as the node stored at a give
% depth and breadth.
% ARG tree : the tree for which the dimensions are required.
% RETURN widths : stores the width at each depth level. 
% RETURN maxDepth : the maximum depth of the tree.
% RETURN nodePositions : stores the nodeIndex present at depth i
% and breadth j.
%
% COPYRIGHT : Andrew J. Moore, 2006
% 
% SEEALSO : treeFindParents, treeFindChildren

% NDLUTIL

maxDepth = 0;
%widths stores the width at each depth level. Since the max depth isn't yet
%known, allocate enough memory to widths for the worst case scenario where
%each node has only 1 child and the depth is equal to the number of nodes.
widths = zeros(length(tree), 1);
%nodePositions 
nodePositions = zeros(length(tree), length(tree));
rootInd = treeFindRoots(tree);
indAlreadySeen = [];
for i = 1:length(rootInd)
  traverseTree(rootInd(i), 1);
end
widths = widths(1:maxDepth, :); %trim off excess rows

  function traverseTree(nodeIndex, depthLevel)
    if ~any(indAlreadySeen == nodeIndex)
      widths(depthLevel) = widths(depthLevel) + 1;
    end
    indAlreadySeen = [indAlreadySeen nodeIndex];

    nodePositions(depthLevel, widths(depthLevel)) = nodeIndex;
    if length(tree(nodeIndex).children) > 0
      for i=1:length(tree(nodeIndex).children)
        traverseTree(tree(nodeIndex).children(i), depthLevel + 1);
      end
    else
      if depthLevel > maxDepth
        maxDepth = depthLevel;
      end
    end
  end

end
