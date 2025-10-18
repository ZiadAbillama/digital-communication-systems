% Test: test_entropy.m
% Purpose:
%   Verify entropy calculation for simple and quantized data.

clc; clear; close all;

% Example 1: simple discrete source
x1 = [1 1 1 2 2 3 3 3 3 3]; % skewed distribution
H1 = entropy(x1);
fprintf('Entropy of x1 = %.4f bits/symbol\n', H1);

% Example 2: uniform distribution
x2 = randi([0 7], [1, 10000]); % 8 equally likely symbols
H2 = entropy(x2);
fprintf('Entropy of x2 = %.4f bits/symbol (expected ~3)\n', H2);

% Example 3: quantized signal
t = 0:0.001:1;
x = sin(2*pi*3*t);
xmin = min(x);
xmax = max(x);
M = 8;
thr = linspace(xmin, xmax, M+1); thr = thr(2:end-1);
lvl = linspace(xmin+(xmax-xmin)/(2*M), xmax-(xmax-xmin)/(2*M), M);
xq = quan(x, thr, lvl);
H3 = entropy(xq);
fprintf('Entropy of quantized signal = %.4f bits/symbol\n', H3);

disp('✅ test_entropy.m executed successfully.');
