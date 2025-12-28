% Purpose:
%   Analyze MSE vs number of quantization levels M using correct uniform quantizer.

t = 0:0.001:1;
x = sin(2*pi*3*t);
x = x(:);

M_values = [2, 4, 8, 16, 32, 64, 128]; %number of representation points
MSE_values = zeros(size(M_values));
MSE_theory = zeros(size(M_values));

xmin = min(x);
xmax = max(x);

for i = 1:length(M_values)
    M = M_values(i);

    % Correct uniform quantizer
    delta = (xmax - xmin) / M;

    thr = xmin + delta * (1:M-1);               % M-1 thresholds
    lvl = xmin + delta*(0:M-1) + delta/2;       % M levels at bin centers

    % Quantize
    xq = quan(x, thr, lvl);

    % MSE
    MSE_values(i) = mean((x - xq).^2);

    % Theoretical MSE
    MSE_theory(i) = delta^2 / 12;
end

%Plot
figure;
semilogy(M_values, MSE_values, 'o-m', 'LineWidth', 1.5); hold on;
semilogy(M_values, MSE_theory, 's--b', 'LineWidth', 1.5);
grid on;

xlabel('Number of Quantization Levels (M)');
ylabel('Mean Squared Error (MSE)');
title('Quantization MSE vs Number of Levels (Uniform Quantizer)');
legend('Measured MSE', 'Theoretical MSE', 'Location', 'southwest');

% Display table
disp(table(M_values.', MSE_values.', MSE_theory.', ...
    'VariableNames', {'M', 'Measured_MSE', 'Theoretical_MSE'}));

disp('test_quan_mse.m executed successfully.');
