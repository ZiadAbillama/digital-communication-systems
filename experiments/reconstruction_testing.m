% RECONSTRUCT_TESTING - Compare reconstructed and original signals
% EECE 442 - Part I: Sampling & Reconstruction

clear; close all; clc;

% -------------------------------------------------------------------------
% Signal Definition
% -------------------------------------------------------------------------
f0 = 100;                % Frequency (Hz)
t = 0:1e-5:0.05;         % Continuous time vector
x = cos(2*pi*f0*t);      % Original signal

fN = 2*f0;               % Nyquist frequency
fs_below = 0.5 * fN;     % Below Nyquist (Aliasing)
fs_above = 2 * fN;       % Above Nyquist (Proper sampling)

% -------------------------------------------------------------------------
% Case 1: Below Nyquist
% -------------------------------------------------------------------------
[tb, xb] = sample(t, x, fs_below);
xr_below = reconstruct(t, xb, fs_below);

figure;
plot(t, x, 'b', 'LineWidth', 1.2); hold on;          % Original in blue
plot(t, xr_below, 'r--', 'LineWidth', 1.2);          % Reconstructed in red dashed
title(sprintf('Reconstruction Below Nyquist (fs = %.1f Hz)', fs_below));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Original Signal', 'Reconstructed Signal');
grid on;
ylim([-1.5 1.5]);                                     % FIXED Y-axis range

% -------------------------------------------------------------------------
% Case 2: Above Nyquist
% -------------------------------------------------------------------------
[ta, xa] = sample(t, x, fs_above);
xr_above = reconstruct(t, xa, fs_above);

figure;
plot(t, x, 'b', 'LineWidth', 1.2); hold on;          % Original in blue
plot(t, xr_above, 'r--', 'LineWidth', 1.2);          % Reconstructed in red dashed
title(sprintf('Reconstruction Above Nyquist (fs = %.1f Hz)', fs_above));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Original Signal', 'Reconstructed Signal');
grid on;
ylim([-1.5 1.5]);                                     % Same Y-axis range

% -------------------------------------------------------------------------
% Error Metrics
% -------------------------------------------------------------------------
err_below = max(abs(x - xr_below));
err_above = max(abs(x - xr_above));

fprintf('Max reconstruction error (below Nyquist): %.3e\n', err_below);
fprintf('Max reconstruction error (above Nyquist): %.3e\n', err_above);

disp('> Below Nyquist: Reconstruction distorted (aliasing).');
disp('> Above Nyquist: Accurate reconstruction (Nyquist satisfied).');
