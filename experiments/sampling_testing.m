% Signal Definition
t = linspace(0, 1, 1000);
xt = cos(2*pi*3*t) + 0.5*cos(2*pi*7*t);

f_max = 7;  % Maximum frequency in the signal

%Under-sampling (fs = 0.5 * f_max = 3.5 Hz)
fs1 = 0.5 * f_max;
[t_sample1, x_sample1] = sample(t, xt, fs1);

figure;
plot(t, xt, 'b', 'LineWidth', 1.5); hold on;
stem(t_sample1, x_sample1, 'r', 'filled');
xlabel('Time (s)'); ylabel('Amplitude');
title('Under-sampling at 0.5 × f_{max} = 3.5 Hz');
legend('Original Signal', 'Sampled Points');
grid on;

%Aliasing Zone (fs = f_max = 7 Hz)
fs2 = f_max;
[t_sample2, x_sample2] = sample(t, xt, fs2);

figure;
plot(t, xt, 'b', 'LineWidth', 1.5); hold on;
stem(t_sample2, x_sample2, 'r', 'filled');
xlabel('Time (s)'); ylabel('Amplitude');
title('Sampling at f_{max} = 7 Hz (Aliasing)');
legend('Original Signal', 'Sampled Points');
grid on;

%Nyquist Rate (fs = 2 * f_max = 14 Hz)
fs3 = 2 * f_max;
[t_sample3, x_sample3] = sample(t, xt, fs3);

figure;
plot(t, xt, 'b', 'LineWidth', 1.5); hold on;
stem(t_sample3, x_sample3, 'r', 'filled');
xlabel('Time (s)'); ylabel('Amplitude');
title('Critical Sampling at 2 × f_{max} = 14 Hz (Nyquist Rate)');
legend('Original Signal', 'Sampled Points');
grid on;

% Oversampling (fs = 4 * f_max = 28 Hz)
fs4 = 4 * f_max;
[t_sample4, x_sample4] = sample(t, xt, fs4);

figure;
plot(t, xt, 'b', 'LineWidth', 1.5); hold on;
stem(t_sample4, x_sample4, 'r', 'filled');
xlabel('Time (s)'); ylabel('Amplitude');
title('Over-sampling at 4 × f_{max} = 28 Hz');
legend('Original Signal', 'Sampled Points');
grid on;
