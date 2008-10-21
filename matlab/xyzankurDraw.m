function handle = xyzankurDraw(joint,handle)

% XYZANKURDRAW
%
% COPYRIGHT : Carl Henrik Ek and Neil Lawrence, 2008

% VISUALISATION


limb{1} = [2 3;3 4;4 5]; %spine
limb{2} = [3 6;6 7;7 8;8 9]; % left-arm
limb{3} = [3 10;10 11;11 12;12 13]; % right-arm
limb{4} = [2 14;14 15;15 16;16 20]; % left-leg
limb{5} = [2 17;17 18;18 19;19 21]; % right-leg


if(nargin<2)
  % draw figure
  k = 1;
  for(i = 1:1:length(limb))
    if(ismember(i,[1 2 4]))
      % left side
      linestyle = '-';
    end
    if(ismember(i,[3 5]))
      % right side
      linestyle = '-.';
    end
    for(j = 1:1:size(limb{i},1))
      handle(k) = line(joint(limb{i}(j,:),1),joint(limb{i}(j,:),2),joint(limb{i}(j,:),3),'LineWidth',3,'LineStyle',linestyle);
      k = k + 1;
    end
  end
else
  % modify figure
  k = 1;
  for(i = 1:1:length(limb))
    for(j = 1:1:size(limb{i},1))
      set(handle(k),'Xdata',joint(limb{i}(j,:),1),'Ydata',joint(limb{i}(j,:),2),'Zdata',joint(limb{i}(j,:),3));
      k = k+1;
    end
  end
end

return;