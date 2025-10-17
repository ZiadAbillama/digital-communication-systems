% Test: test_lloyd_max.m
% Purpose: Compare uniform vs Lloyd–Max quantization for a non-uniform signal.

clc; clear; close all;

% Generate a non-uniform signal (Laplace or exponential)
x = randn(1, 10000); % Gaussian-distributed samples

M = 8; % number of quantization levels

% Lloyd–Max quantizer
[xq_opt, thr_opt, lvl_opt] = lloyd_max(x, M);

% Uniform quantizer for comparison
xmin = min(x);
xmax = max(x);
thr_uni = linspace(xmin, xmax, M+1);
thr_uni = thr_uni(2:end-1);
lvl_uni = linspace(xmin + (xmax-xmin)/(2*M), xmax - (xmax-xmin)/(2*M), M);
xq_uni = quan(x, thr_uni, lvl_uni);

% Compute MSE
MSE_uni = mean((x - xq_uni).^2);
MSE_opt = mean((x - xq_opt).^2);

% Display comparison
fprintf('Uniform Quantizer MSE:     %.6f\n', MSE_uni);
fprintf('Lloyd–Max Quantizer MSE:  %.6f\n', MSE_opt);

% Plot histograms for visualization
figure;
subplot(2,1,1);
histogram(x, 100, 'Normalization', 'pdf'); hold on;
stem(lvl_uni, zeros(size(lvl_uni)), 'r', 'LineWidth', 1.5);
title('Uniform Quantization Levels');
legend('Input PDF', 'Levels');

subplot(2,1,2);
histogram(x, 100, 'Normalization', 'pdf'); hold on;
stem(lvl_opt, zeros(size(lvl_opt)), 'm', 'LineWidth', 1.5);
title('Lloyd–Max Quantization Levels');
legend('Input PDF', 'Optimized Levels');

disp('✅ test_lloyd_max.m executed successfully');
