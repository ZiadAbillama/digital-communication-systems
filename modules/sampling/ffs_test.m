% ffs_test.m  –  Testing the Finite Fourier Series approximation
clc; clear; close all;

% Define time-limited signal
T = 1;                              % fundamental period
t = linspace(-T/2, T/2, 1000);
xt = cos(2*pi*3*t) + 0.5*cos(2*pi*7*t);   % example message chunk

% Test different numbers of harmonics n
n_values = [3, 7, 15, 30];
figure;
for i = 1:length(n_values)
    [xhat, ~] = ffs(xt, t, n_values(i), T);
    subplot(2,2,i);
    plot(t, xt, 'b', t, xhat, 'r--');
    title(['n = ', num2str(n_values(i))]);
    xlabel('t'); ylabel('x(t)');
    legend('Original','Fourier approx');
    grid on;
end
sgtitle('Finite Fourier Series Approximations for different n');

% ---------------------------------------------------------------
% error analysis  E(n,T) vs n
E = zeros(size(n_values));
for i = 1:length(n_values)
    [xhat, ~] = ffs(xt, t, n_values(i), T);
    E(i) = trapz(t, abs(xt - xhat).^2);      % numerical integral
end
figure;
plot(n_values, E, '-o');
xlabel('n (number of harmonics)');
ylabel('Error energy E(n,T)');
title('Approximation error vs n');
grid on;

% ---------------------------------------------------------------
% effect of varying T (keep n large)
T_values = [0.5, 1, 2];
n_large = 30;
E_T = zeros(size(T_values));
for i = 1:length(T_values)
    [xhat, ~] = ffs(xt, t, n_large, T_values(i));
    E_T(i) = trapz(t, abs(xt - xhat).^2);
end
figure;
plot(T_values, E_T, '-s');
xlabel('T');
ylabel('Error energy E(n,T)');
title('Approximation error vs T (for large n)');
grid on;
