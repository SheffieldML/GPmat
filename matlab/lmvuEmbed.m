function X = lmvuEmbed(Y,dims,k,nr_landmark)

% LMVUEMBED Embed data set with landmark MVU
% FORMAT
% DESC Embed data set with landmark version of Weinberg et al.'s
% maximum variance unfolding algorithm.
% ARG Y : Data
% ARG dims : Dimensionality of Embedding (default = 2)
% ARG k : Number of Neighbours in Proximity Graph (default = 7)
% ARG nr_landmark : Number of landmark Points
% RETURN X : embedding
%
% SEEALSO : ppcaEmbed, lleEmbed, mvuEmbed
%
% COPYRIGHT : Carl Henrik Ek, Neil D. Lawrence, 2007

% MLTOOLS

if(nargin<4)
  nr_landmark = 30;
  if(nargin<3)
    k = 7;
    if(nargin<2)
      dims = 2;
      if(nargin<1)
	error('To Few Arguments');
      end
    end
  end
end
  
if(any(any(isnan(Y))))
  error('Cannot Initialise GPLVM using lmvu when missing data is present.');
end

X = lmvu(distance(Y'),nr_landmark,k);

X = X(1:1:dims,:)';

return
