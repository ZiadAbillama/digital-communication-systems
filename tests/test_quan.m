% Test: test_quan.m
% Purpose: Verify basic two-level and uniform multi-level quantization.

clc; clear; close all;

% Generate a test signal
t = 0:0.001:1;
x = sin(2*pi*3*t); % Simple sine wave

% Two-level quantization
thr = 0;           % threshold at zero
lvl = [-1 1];      % two quantization levels

xq_2 = quan(x, thr, lvl);

% Uniform multi-level quantization (M = 8)
M = 8;
xmin = min(x);
xmax = max(x);
thr_uniform = linspace(xmin, xmax, M+1);
thr_uniform = thr_uniform(2:end-1); % remove first and last
lvl_uniform = linspace(xmin + (xmax-xmin)/(2*M), xmax - (xmax-xmin)/(2*M), M);

xq_8 = quan(x, thr_uniform, lvl_uniform);

% Plot comparison
figure;
subplot(2,1,1);
plot(t, x, 'b', 'LineWidth', 1.2); hold on;
stairs(t, xq_2, 'r', 'LineWidth', 1.2);
title('Two-Level Quantization');
legend('Original', 'Quantized');
grid on;

subplot(2,1,2);
plot(t, x, 'b', 'LineWidth', 1.2); hold on;
stairs(t, xq_8, 'm', 'LineWidth', 1.2);
title('Uniform 8-Level Quantization');
legend('Original', 'Quantized');
grid on;

disp('✅ test_quan.m executed successfully');
