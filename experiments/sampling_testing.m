% SAMPLE_TESTING - Visualizes sampling below and above Nyquist
% EECE 442 - Part I

clear; close all; clc;

% -------------------------------------------------------------------------
% Signal Definition
% -------------------------------------------------------------------------
f0 = 100;                % Signal frequency (Hz)
t = 0:1e-5:0.05;         % "Continuous" reference time
x = cos(2*pi*f0*t);      % Original signal

fN = 2*f0;               % Nyquist frequency
fs_below = 0.5 * fN;     % Below Nyquist
fs_above = 2 * fN;       % Above Nyquist

% -------------------------------------------------------------------------
% Case 1: Below Nyquist
% -------------------------------------------------------------------------
[tb, xb] = sample(t, x, fs_below);

figure;
plot(t, x, 'b', 'LineWidth', 1.2); hold on;        % Original signal in blue
stem(tb, xb, 'r', 'filled');                       % Samples in red
title(sprintf('Sampling Below Nyquist (fs = %.1f Hz)', fs_below));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Original Signal', 'Samples');
grid on;

% -------------------------------------------------------------------------
% Case 2: Above Nyquist
% -------------------------------------------------------------------------
[ta, xa] = sample(t, x, fs_above);

figure;
plot(t, x, 'b', 'LineWidth', 1.2); hold on;        % Original signal in blue
stem(ta, xa, 'r', 'filled');                       % Samples in red
title(sprintf('Sampling Above Nyquist (fs = %.1f Hz)', fs_above));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Original Signal', 'Samples');
grid on;

disp('> Below Nyquist: Aliasing visible (samples too sparse)');
disp('> Above Nyquist: Proper sampling (samples follow signal)');
