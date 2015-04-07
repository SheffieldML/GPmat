t = -3:.1:3;

[X, Y] = meshgrid(t, t);

c = 0.95;

C = toeplitz([1, c]);
L = inv(C);

Z = exp(- (L(1, 1) * X .^ 2 + L(2, 2) * Y .^ 2 + 2 * L(1, 2) * X .* Y));

figure(1);

subplot(4, 4, [2:4, 6:8, 10:12]);
contour(X, Y, Z);
set(gca, 'XTick', []);
set(gca, 'YTick', []);

subplot(4, 4, [1 5 9]);
plot(1 ./ sqrt(2*pi) * exp(- .5 * t.^2), t);
set(gca, 'XTick', []);
set(gca, 'YTick', []);

subplot(4, 4, 14:16);
plot(t, 1 ./ sqrt(2*pi) * exp(- .5 * t.^2));
set(gca, 'XTick', []);
set(gca, 'YTick', []);

figure(2);

subplot(4, 4, [2:4, 6:8, 10:12]);
contour(X, Y, Z);
hold on
plot([-3, 3], [1.5, 1.5])
set(gca, 'XTick', []);
set(gca, 'YTick', []);
hold off

subplot(4, 4, [1 5 9]);
plot(1 ./ sqrt(2*pi) * exp(- .5 * t.^2), t, '--');
hold on
plot(1, 1.5, '*');
plot([0 1], [1.5, 1.5]);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
hold off
    
subplot(4, 4, 14:16);
plot(t, 1 ./ sqrt(2*pi) * exp(- .5 * t.^2), '--');
hold on
plot(t, 1 ./ sqrt(2*pi*(1-c.^2)) * exp(- .5 ./ (1 - c.^2) * (t - c * 1.5).^2));
set(gca, 'XTick', []);
set(gca, 'YTick', []);
hold off

