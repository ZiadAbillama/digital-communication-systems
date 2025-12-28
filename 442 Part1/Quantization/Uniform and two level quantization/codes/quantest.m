% Purpose:
%   Verify two-level and uniform multi-level quantization using the correct formulas.

t = 0:0.001:1;
x = sin(2*pi*3*t);   % Simple test signal
x = x(:);

% Two-Level Quantization 
thr2 = 0;                 % Threshold at zero
lvl2 = [-1, 1];           % Two representation points
xq2 = quan(x, thr2, lvl2);

% Uniform Multi-Level Quantization
M = 32;                    % Number of levels
xmin = min(x);
xmax = max(x);

% Correct uniform quantizer:
delta = (xmax - xmin) / M;

thr_uni = xmin + delta * (1:M-1);          % thresholds
lvl_uni = xmin + delta*(0:M-1) + delta/2;  % representation points

xq_uni = quan(x, thr_uni, lvl_uni);

%Plot
figure;
subplot(2,1,1);
plot(t, x, 'b', 'LineWidth', 1.2); hold on;
stairs(t, xq2, 'r', 'LineWidth', 1.2);
title('Two-Level Quantization');
legend('Original', 'Quantized');
grid on;

subplot(2,1,2);
plot(t, x, 'b', 'LineWidth', 1.2); hold on;
stairs(t, xq_uni, 'm', 'LineWidth', 1.2);
title(sprintf('Uniform %d-Level Quantization (Corrected)', M));
legend('Original', 'Quantized');
grid on;

disp('test_quan.m executed successfully.');
