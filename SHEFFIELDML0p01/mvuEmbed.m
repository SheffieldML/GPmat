function X = mvuEmbed(Y,dims,k)

% MVUEMBED Embed data set with MVU.
%
%	Description:
%
%	X = MVUEMBED(Y, DIMS, K) Embed a given data set with Weinberger et
%	al.'s MVU algorithm.
%	 Returns:
%	  X - embedding
%	 Arguments:
%	  Y - Data
%	  DIMS - Dimensionality of Embedding (default = 2)
%	  K - Number of Neighbours in Proximity Graph (default = 7)
%	
%
%	See also
%	PPCAEMBED, LLEEMBED, LMVUEMBED


%	Copyright (c) Neil D. Lawrence, 2007 Carl Henrik Ek


if(nargin<3)
  k = 7;
  if(nargin<2)
    dims = 2;
    if(nargin<1)
      error('To Few Arguments');
    end
  end
end

if(any(any(isnan(Y))))
  error('Cannot run MVU when missing data is present.');
end

X = mvu(distance(Y'),k);

X = X(1:1:dims,:)';


function D = distance(Y)
  
  D = sqrt(dist2(Y', Y'));
return