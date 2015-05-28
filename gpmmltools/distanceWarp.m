function [w, Dist,D,k]=distanceWarp(d)

% DISTANCEWARP Dynamic Time Warping Algorithm
%
% DESC returns the path through a pair of dynamic time warped signals.
% ARG d : distance between the two signals.
% RETURN w : is the optimal path
% RETURN Dist : is unnormalized distance between the signals.
% RETURN D : is the accumulated distance matrix
% RETURN k : is the normalizing factor

% MLTOOLS
  
  [N, M] = size(d);

  D=zeros(size(d));
  D(1,1)=d(1,1);
  
  for n=2:N
    D(n,1)=d(n,1)+D(n-1,1);
  end
  for m=2:M
    D(1,m)=d(1,m)+D(1,m-1);
  end
  for n=2:N
    for m=2:M
      D(n,m)=d(n,m)+min([D(n-1,m),D(n-1,m-1),D(n,m-1)]);
    end
  end
  
  Dist=D(N,M);
  n=N;
  m=M;
  k=1;
  w=[];
  w(1,:)=[N,M];
  while ((n+m)~=2)
    if (n-1)==0
      m=m-1;
    elseif (m-1)==0
      n=n-1;
    else 
      [values,number]=min([D(n-1,m),D(n,m-1),D(n-1,m-1)]);
      switch number
       case 1
        n=n-1;
       case 2
        m=m-1;
       case 3
        n=n-1;
        m=m-1;
      end
    end
    k=k+1;
    w=cat(1,w,[n,m]);
  end
end
