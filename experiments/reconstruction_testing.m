% Original signal
t = linspace(0, 1, 1000);
xt = cos(2*pi*3*t) + 0.5*cos(2*pi*7*t);
f_max = 7;  % Maximum frequency component

t_rec = t; % reconstruction grid

%% Under-sampling: fs = 0.5 × fmax = 3.5 Hz
fs_under = 0.5 * f_max;
[t_sample_under, x_sample_under] = sample(t, xt, fs_under);
[~, x_rec_under] = reconstruct(t_sample_under, x_sample_under, t_rec);

figure;
plot(t, xt, 'b'); hold on;
plot(t_rec, x_rec_under, 'r--');
legend('Original Signal', 'Reconstructed Signal');
xlabel('Time (s)'); ylabel('Amplitude');
title('Reconstruction after Under-sampling (0.5 × f_{max} = 3.5 Hz)');
grid on;

%% Aliasing Case: fs = fmax = 7 Hz
fs_alias = f_max;
[t_sample_alias, x_sample_alias] = sample(t, xt, fs_alias);
[~, x_rec_alias] = reconstruct(t_sample_alias, x_sample_alias, t_rec);

figure;
plot(t, xt, 'b'); hold on;
plot(t_rec, x_rec_alias, 'r--');
legend('Original Signal', 'Reconstructed Signal');
xlabel('Time (s)'); ylabel('Amplitude');
title('Reconstruction after Sampling at f_{max} = 7 Hz (Aliasing)');
grid on;

%% Critical Sampling: fs = 2 × fmax = 14 Hz
fs_nyquist = 2 * f_max;
[t_sample_nyq, x_sample_nyq] = sample(t, xt, fs_nyquist);
[~, x_rec_nyq] = reconstruct(t_sample_nyq, x_sample_nyq, t_rec);

figure;
plot(t, xt, 'b'); hold on;
plot(t_rec, x_rec_nyq, 'r--');
legend('Original Signal', 'Reconstructed Signal');
xlabel('Time (s)'); ylabel('Amplitude');
title('Reconstruction at Critical Sampling (2 × f_{max} = 14 Hz)');
grid on;

%% Over-sampling: fs = 4 × fmax = 28 Hz
fs_over = 4 * f_max;
[t_sample_over, x_sample_over] = sample(t, xt, fs_over);
[~, x_rec_over] = reconstruct(t_sample_over, x_sample_over, t_rec);

figure;
plot(t, xt, 'b'); hold on;
plot(t_rec, x_rec_over, 'r--');
legend('Original Signal', 'Reconstructed Signal');
xlabel('Time (s)'); ylabel('Amplitude');
title('Reconstruction after Over-sampling (4 × f_{max} = 28 Hz)');
grid on;

%% Extra: Gaussian Pulse
t = linspace(0, 1, 1000);
xt = exp(-10*(t-0.5).^2);
fs_gauss = 30;
[t_sample_gauss, x_sample_gauss] = sample(t, xt, fs_gauss);
[~, x_rec_gauss] = reconstruct(t_sample_gauss, x_sample_gauss, t);

figure;
plot(t, xt, 'b'); hold on;
plot(t, x_rec_gauss, 'r--');
legend('Original Signal', 'Reconstructed Signal');
xlabel('Time (s)'); ylabel('Amplitude');
title('Reconstruction of Gaussian Pulse');
grid on;
