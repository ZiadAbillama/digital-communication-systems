% Test: test_reconstruct.m
% Purpose: Verify that reconstruct.m correctly rebuilds the signal.

clc; clear; close all;

% Original signal
t = 0:0.0001:1;
xt = cos(2*pi*10*t);

% Sampling
fs = 20; % Hz
[t_sample, x_sample] = sample(t, xt, fs);

% Reconstruction
xrcon = reconstruct(t, x_sample, fs);

% Plot comparison
figure;
plot(t, xt, 'b', 'LineWidth', 1.2); hold on;
plot(t, xrcon, 'r--', 'LineWidth', 1.2);
legend('Original', 'Reconstructed');
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal Reconstruction Using sinc Interpolation');
grid on;

disp('✅ test_reconstruct.m executed successfully');
