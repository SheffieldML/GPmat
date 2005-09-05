function connection = bvhConnectionMatrix(bvhStruct);

% BVHCONNECTIONMATRIX Compute the connection matrix for the structure.

% MOCAP

connection = zeros(length(bvhStruct));
for i = 1:length(bvhStruct);
  for j = 1:length(bvhStruct(i).children)    
    connection(i, bvhStruct(i).children(j)) = 1;
  end
end

