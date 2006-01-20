function tree = treeSwapNode(tree, i, j);

% TREESWAPNODE Swap two nodes in the tree structure array.

% MOCAP

storeNodeI = tree(i);
storeNodeJ = tree(j);
tree(j) = storeNode(i);
tree(i) = storeNode(j);
for k = 1:length(tree)
  tree(k).children(find(tree(k).children)==i) = -1;
  tree(k).children(find(tree(k).children)==j) = i;
  tree(k).children(find(tree(k).children)==-1) = j;
  tree(k).parent(find(tree(k).parent)==i) = -1;
  tree(k).parent(find(tree(k).parent)==j) = i;
  tree(k).parent(find(tree(k).parent)==-1) = j;
end
