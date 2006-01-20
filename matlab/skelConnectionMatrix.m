function connection = skelConnectionMatrix(skel);

% SKELCONNECTIONMATRIX Compute the connection matrix for the structure.

% MOCAP

connection = zeros(length(skel.tree));
for i = 1:length(skel.tree);
  for j = 1:length(skel.tree(i).children)    
    connection(i, skel.tree(i).children(j)) = 1;
  end
end

