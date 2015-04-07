function plot_gp(rho, style),

switch style,
 case 1,
  subplot(2, 2, 1);
  set(gca, 'FontSize', 16);
  errorbar([0, 1], [0, 0], [2, 2], 'x');
  axis([-1 2 -2.5 2.5]);
  %set(gca, 'XTick', []);
  %set(gca, 'YTick', []);

  subplot(2, 2, 2);
  set(gca, 'FontSize', 16);
  obs = 1.5;
  m = rho * obs;
  v = 1 - rho.^2;
  errorbar([0, 1], [obs, m], [0, 2*sqrt(v)], 'x');
  axis([-1 2 -2.5 2.5]);
  %set(gca, 'XTick', []);
  %set(gca, 'YTick', []);
 case 2,
  subplot(2, 2, 1);
  set(gca, 'FontSize', 16);

  plot([-1, 2], [0, 0]);
  hold on;
  plot([-1, 2], [2, 2], '--');
  plot([-1, 2], [-2, -2], '--');
  errorbar([0, 1], [0, 0], [2, 2], 'x');
  hold off;
  axis([-1 2 -2.5 2.5]);
  %set(gca, 'XTick', []);
  %set(gca, 'YTick', []);

  subplot(2, 2, 2);
  set(gca, 'FontSize', 16);
  obs = 1.5;
  t = -1:.01:2;
  l2 = -1 / log(rho);
  
  m1 = rho * obs;
  v1 = 1 - rho.^2;
  
  k = exp(-t.^2/l2);
  m = k * obs;
  v = 1 - k.^2;
  plot(t, m);
  hold on;
  plot(t, m + 2*sqrt(v), '--');
  plot(t, m - 2*sqrt(v), '--');
  errorbar([0, 1], [obs, m1], [0, 2*sqrt(v1)], 'x');
  hold off;
  axis([-1 2 -2.5 2.5]);
  %set(gca, 'XTick', []);
  %set(gca, 'YTick', []);
end
