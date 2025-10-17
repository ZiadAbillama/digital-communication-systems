% Test: test_ffs_error.m
% Purpose:
%   Compute and plot the Fourier Series approximation error E(n, T)
%   both as a function of n (number of harmonics) and T (period).

clc; clear; close all;

% Define time and original signal
t = linspace(-1, 1, 2000);
xt = square(2*pi*2*t); % 2 Hz square wave signal

% Case 1: Fixed T, vary n
T_fixed = 2;                         % period (seconds)
n_values = [1, 10, 20, 40, 80];
E_n = zeros(size(n_values));

for i = 1:length(n_values)
    n = n_values(i);
    [xhat, ~] = ffs(xt, t, n, T_fixed);
    E_n(i) = trapz(t, abs(xt - xhat).^2);
end

% Case 2: Fixed n, vary T
n_fixed = 20;                        % number of harmonics
T_values = [0.5, 1, 2, 4, 8];
E_T = zeros(size(T_values));

for i = 1:length(T_values)
    T = T_values(i);
    [xhat, ~] = ffs(xt, t, n_fixed, T);
    E_T(i) = trapz(t, abs(xt - xhat).^2);
end

% Plot both results
figure('Name','Fourier Series Error Analysis','NumberTitle','off');

subplot(1,2,1);
semilogy(n_values, E_n, 'or-', 'LineWidth', 1.5);
grid on;
xlabel('Number of Harmonics (n)');
ylabel('Error E(n, T)');
title(sprintf('Error vs n (T = %.1f s)', T_fixed));
legend('E(n, T)', 'Location', 'northeast');

subplot(1,2,2);
semilogy(T_values, E_T, 'ob-', 'LineWidth', 1.5);
grid on;
xlabel('Approximation Period (T) [s]');
ylabel('Error E(n, T)');
title(sprintf('Error vs T (n = %d)', n_fixed));
legend('E(n, T)', 'Location', 'northeast');

sgtitle('Finite Fourier Series Approximation Error Analysis');

% Display numeric results
disp('--- Error vs n (T fixed) ---');
disp(table(n_values.', E_n.', 'VariableNames', {'n', 'Error_E'}));

disp('--- Error vs T (n fixed) ---');
disp(table(T_values.', E_T.', 'VariableNames', {'T', 'Error_E'}));

disp('✅ test_ffs_error.m executed successfully');
