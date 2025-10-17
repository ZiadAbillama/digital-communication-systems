% Test: test_ffs.m
% Purpose: Verify finite Fourier series approximation.

clc; clear; close all;

% Define time and signal
t = linspace(-1, 1, 2000);
xt = square(2*pi*2*t); % periodic signal with discontinuities

% Parameters
n = 10; % number of harmonics
T = 2;  % fundamental period

[xhat, ck] = ffs(xt, t, n, T);

% Plot the original and approximated signals
figure;
plot(t, xt, 'b', 'LineWidth', 1.2); hold on;
plot(t, xhat, 'r--', 'LineWidth', 1.2);
legend('Original', sprintf('FFS Approximation (n=%d)', n));
xlabel('Time (s)');
ylabel('Amplitude');
title('Finite Fourier Series Approximation');
grid on;

disp('✅ test_ffs.m executed successfully');
