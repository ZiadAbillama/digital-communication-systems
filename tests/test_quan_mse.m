% Test: test_quan_mse.m
% Purpose:
%   Analyze how Mean Squared Error (MSE) changes with number of quantization levels M.

clc; clear; close all;

% Generate a test signal
t = 0:0.001:1;
x = sin(2*pi*3*t); % 3 Hz sine wave
x = x(:);

% Range of quantization levels to test
M_values = [2, 4, 8, 16, 32, 64, 128];
MSE_values = zeros(size(M_values));

% Loop over each quantization level
for i = 1:length(M_values)
    M = M_values(i);
    
    % Compute thresholds and levels for uniform quantization
    xmin = min(x);
    xmax = max(x);
    thr = linspace(xmin, xmax, M+1);
    thr = thr(2:end-1); % remove first and last points
    lvl = linspace(xmin + (xmax-xmin)/(2*M), xmax - (xmax-xmin)/(2*M), M);
    
    % Quantize the signal
    xq = quan(x, thr, lvl);
    
    % Compute mean squared error
    MSE_values(i) = mean((x - xq).^2);
end

% Plot MSE vs number of levels
figure;
semilogy(M_values, MSE_values, 'o-m', 'LineWidth', 1.5);
grid on;
xlabel('Number of Quantization Levels (M)');
ylabel('Mean Squared Error (MSE)');
title('Quantization Error vs. Number of Levels');
legend('MSE', 'Location', 'northeast');

% Print numerical results
disp(table(M_values.', MSE_values.', 'VariableNames', {'Levels_M', 'MSE'}));

disp('✅ test_quan_mse.m executed successfully');
