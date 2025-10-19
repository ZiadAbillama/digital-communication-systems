% RECONSTRUCT_TESTING - Compare reconstructed and original signals

% Signal Definition
f0 = 100;                % Frequency (Hz)
t = 0:1e-5:0.05;         % Continuous time vector
x = cos(2*pi*f0*t);      % Original signal

fN = 2*f0;               % Nyquist frequency
fs_below = 0.5 * fN;     % Below Nyquist (Aliasing)
fs_above = 2 * fN;       % Above Nyquist (Proper sampling)

% Case 1: Below Nyquist
[tb, xb] = sample(t, x, fs_below);
xr_below = reconstruct(t, xb, fs_below);

figure;
plot(t, x, 'b', 'LineWidth', 1.2); hold on;          % Original 
plot(t, xr_below, 'r--', 'LineWidth', 1.2);          % Reconstructed
title(sprintf('Reconstruction Below Nyquist (fs = %.1f Hz)', fs_below));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Original Signal', 'Reconstructed Signal');
grid on;
ylim([-1.5 1.5]);                                    


% Case 2: Above Nyquist
[ta, xa] = sample(t, x, fs_above);
xr_above = reconstruct(t, xa, fs_above);

figure;
plot(t, x, 'b', 'LineWidth', 1.2); hold on;          % Original
plot(t, xr_above, 'r--', 'LineWidth', 1.2);          % Reconstructed
title(sprintf('Reconstruction Above Nyquist (fs = %.1f Hz)', fs_above));
xlabel('Time (s)'); ylabel('Amplitude');
legend('Original Signal', 'Reconstructed Signal');
grid on;
ylim([-1.5 1.5]);  

% MSE
mse = @(a,b) mean((a(:) - b(:)).^2);

mse_below = mse(x, xr_below);
mse_above = mse(x, xr_above);

fprintf('MSE (below Nyquist): %.3e\n', mse_below);
fprintf('MSE (above Nyquist): %.3e\n', mse_above);
