% test_ffs.m
n = 50;       % number of harmonics on each side
T = 2;        % chosen approximation period (seconds)

% Define time over exactly one analysis period and the signal on it
t  = linspace(-T/2, T/2, 2000);
xt = sign(sin(2*pi*2*t));   % 2 Hz square-like signal (true period = 0.5 s)

% FFS approximation (integrates over [-T/2, T/2] internally)
[xhat, ck] = ffs(xt, t, n, T);

% Plot
figure;
plot(t, xt, 'b', 'LineWidth', 1.2); hold on;
plot(t, xhat, 'r--', 'LineWidth', 1.2);
legend('Original', sprintf('FFS Approximation (n=%d)', n), 'Location','best');
xlabel('Time (s)'); ylabel('Amplitude');
title('Finite Fourier Series Approximation (single-period domain)');
grid on;

disp('test_ffs.m executed successfully');
