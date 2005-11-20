function zeroAxes(axesHandle, tickRatio, fontSize, fontName)

% ZEROAXES A function to move the axes crossing point to the origin.

% NDLUTIL

axis off
xlim = get(axesHandle, 'xlim');
ylim = get(axesHandle, 'ylim');
xticklength = (xlim(2) - xlim(1))*tickRatio;
yticklength = (ylim(2) - ylim(1))*tickRatio;

xtick = get(axesHandle, 'xtick');
ytick = get(axesHandle, 'ytick');
xaxYpos = 0;
if ylim(1) > 0
  xaxYpos = ylim(1);
elseif ylim(2) < 0
  xaxYpos = ylim(2);
end
line([xlim(1) xlim(2)], [xaxYpos xaxYpos])
for i = 1:length(xtick)
  if xtick(i) == 0
    if any(ytick == 0)
      continue
    end
  end
  line([xtick(i) xtick(i)], xaxYpos+[0 -yticklength]);
  textHandle = text(xtick(i), xaxYpos-4*yticklength, num2str(xtick(i)));
  set(textHandle, 'fontsize', fontSize, 'fontname', fontName, ...
		  'horizontalalignment', 'center');
end

yaxXpos = 0;
if xlim(1) > 0
  yaxXpos = xlim(1);
  set(axesHandle, 'xlim', [xlim(1) - xticklength xlim(2)]);
elseif xlim(2) < 0
  yaxXpos = xlim(2);
  set(axesHandle, 'xlim', [xlim(1) + xticklength xlim(2)]);
end
line([yaxXpos  yaxXpos], [ylim(1) ylim(2)])
for i = 1:length(ytick)
  if ytick(i) == 0
    if any(xtick ==  0)
      origin = ovalCreate([0 0], xticklength, yticklength);
      origin.selected = 0;
      origin = ovalDraw(origin);
      continue
    end
  end
  line(yaxXpos-[0 xticklength], [ytick(i) ytick(i)]);
  textHandle = text(yaxXpos - 1.1*xticklength, ytick(i),  num2str(ytick(i)));
  set(textHandle, 'fontsize', fontSize, 'fontname', fontName, ...
		  'horizontalalignment', 'right');

end




