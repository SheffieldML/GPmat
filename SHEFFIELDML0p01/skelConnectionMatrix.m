function connection = skelConnectionMatrix(skel);

% SKELCONNECTIONMATRIX Compute the connection matrix for the structure.
%
%	Description:
%
%	CONNECTION = SKELCONNECTIONMATRIX(SKEL) computes the connection
%	matrix for the structure. Returns a matrix which has zeros at all
%	entries except those that are connected in the skeleton.
%	 Returns:
%	  CONNECTION - connectivity matrix.
%	 Arguments:
%	  SKEL - the skeleton for which the connectivity is required.
%	
%
%	See also
%	SKELVISUALISE, SKELMODIFY


%	Copyright (c) 2006 Neil D. Lawrence


connection = zeros(length(skel.tree));
for i = 1:length(skel.tree);
  for j = 1:length(skel.tree(i).children)    
    connection(i, skel.tree(i).children(j)) = 1;
  end
end

