% Test: test_sample.m
% Purpose: Verify that sample.m performs correct signal sampling.

clc; clear; close all;

% Generate a simple test signal (10 Hz cosine)
t = 0:0.001:1; 
xt = cos(2*pi*10*t);

% Sampling rate
fs = 20; % Hz (above Nyquist)
[t_sample, x_sample] = sample(t, xt, fs);

% Plot the results
figure;
plot(t, xt, 'b', 'LineWidth', 1.2); hold on;
stem(t_sample, x_sample, 'r', 'filled');
legend('Original', 'Sampled');
xlabel('Time (s)');
ylabel('Amplitude');
title('Sampling Test: sample.m');
grid on;

disp('✅ test_sample.m executed successfully');
