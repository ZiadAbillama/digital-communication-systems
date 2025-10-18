% Test: test_lloyd_comparison.m
% Purpose:
%   Compare Uniform and Lloyd–Max quantization on the same input signal.
%   (Fulfills requirements f, g, and h from the project specification.)

clc; clear; close all;

%% Step 1: Generate or load your empirical input signal
% For now we can use the reconstructed signal from the sampling module,
% or a representative signal (e.g., a sine wave). Later, this should be
% replaced by your actual sampled/reconstructed signal.

fs = 1000;              % Sampling frequency
t = 0:1/fs:1;           % 1 second duration
x = sin(2*pi*3*t) + 0.3*sin(2*pi*7*t);  % Example composite waveform
x = x(:);               % column vector

%% Step 2: Define quantization parameters
M = 8;  % number of quantization levels to compare
fprintf('Comparing Uniform vs Lloyd–Max quantization with M = %d levels...\n', M);

%% Step 3: Uniform Quantization
xmin = min(x);
xmax = max(x);
thr_uni = linspace(xmin, xmax, M+1);
thr_uni = thr_uni(2:end-1);
lvl_uni = linspace(xmin + (xmax-xmin)/(2*M), xmax - (xmax-xmin)/(2*M), M);
xq_uni = quan(x, thr_uni, lvl_uni);

%% Step 4: Lloyd–Max Quantization (Near-optimal)
[xq_opt, thr_opt, lvl_opt] = lloyd_max(x, M);

%% Step 5: Compute MSEs
MSE_uni = mean((x - xq_uni).^2);
MSE_opt = mean((x - xq_opt).^2);

fprintf('Uniform Quantizer MSE:     %.6e\n', MSE_uni);
fprintf('Lloyd–Max Quantizer MSE:  %.6e\n', MSE_opt);

%% Step 6: Plot reconstructed signals
figure('Name', 'Quantization Comparison', 'NumberTitle', 'off');

subplot(2,1,1);
plot(t, x, 'k', 'LineWidth', 1.2); hold on;
stairs(t, xq_uni, 'r', 'LineWidth', 1);
title(sprintf('Uniform Quantization (M = %d) | MSE = %.2e', M, MSE_uni));
xlabel('Time [s]'); ylabel('Amplitude');
legend('Original', 'Uniform Quantized');
grid on;

subplot(2,1,2);
plot(t, x, 'k', 'LineWidth', 1.2); hold on;
stairs(t, xq_opt, 'b', 'LineWidth', 1);
title(sprintf('Lloyd–Max Quantization (M = %d) | MSE = %.2e', M, MSE_opt));
xlabel('Time [s]'); ylabel('Amplitude');
legend('Original', 'Lloyd–Max Quantized');
grid on;

%% Step 7: Plot histograms of levels
figure('Name', 'Quantization Levels Distribution', 'NumberTitle', 'off');
histogram(x, 100, 'Normalization', 'pdf', 'FaceAlpha', 0.4); hold on;
stem(lvl_uni, zeros(size(lvl_uni)), 'r', 'LineWidth', 1.5);
stem(lvl_opt, zeros(size(lvl_opt)), 'b', 'LineWidth', 1.5);
title('Quantization Level Placement');
xlabel('Signal Amplitude'); ylabel('PDF');
legend('Signal PDF', 'Uniform Levels', 'Lloyd–Max Levels');
grid on;

disp('✅ test_lloyd_comparison.m executed successfully.');
