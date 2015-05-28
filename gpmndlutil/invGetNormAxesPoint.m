function point = invGetNormAxesPoint(axesPoint, axesHandle)

% INVGETNORMAXESPOINT Take a point on a plot and return a point within the figure.

% SHEFFIELDML

position = get(axesHandle, 'Position');
xLim = get(axesHandle, 'XLim');
yLim = get(axesHandle, 'YLim');
xSpan = xLim(2) - xLim(1);
ySpan = yLim(2) - yLim(1);

x = axesPoint(1);
y = axesPoint(2);

x = (x - xLim(1))/xSpan;
y = (y - yLim(1))/ySpan;

point(1) = x*position(3) + position(1);
point(2) = y*position(4) + position(2);

