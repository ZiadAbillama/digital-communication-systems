% QPSK over ISI channel with Viterbi equalizer
% Channel: y_i = 2*x_i + x_{i-1} + n_i

clear; close all; clc;

N_bits   = 2e5;
EbN0_dB  = 0:2:14;
N_SNR    = length(EbN0_dB);

A  = 1;          % QPSK amplitude -> Es = 1
Es = 1;
Eb = Es/2;       % 2 bits per QPSK symbol

h = [2 1];                      % ISI taps
h_energy = sum(abs(h).^2);      % ||h||^2 = 5

BER_AWGN    = zeros(1, N_SNR);
BER_ISI_VIT = zeros(1, N_SNR);

% QPSK constellation (same mapping as qpsk_mod / qpsk_demod)
bits_list = [1 1; 1 0; 0 0; 0 1];
S = zeros(1,4);
for i = 1:4
    S(i) = qpsk_mod(bits_list(i,:), A);
end

fprintf('=== QPSK Constellation ===\n');
for i = 1:4
    fprintf('  Symbol %d: [%d %d] -> (%.3f, %.3f)\n', ...
        i, bits_list(i,1), bits_list(i,2), real(S(i)), imag(S(i)));
end
fprintf('\n');

for idx = 1:N_SNR

    EbN0_lin = 10^(EbN0_dB(idx)/10);
    EsN0_lin = 2*EbN0_lin;          % Es/N0 for QPSK
    N0       = Es / EsN0_lin;

    % AWGN noise std
    sigma_awgn = sqrt(N0/2);

    % For ISI channel, scale noise so that effective SNR matches AWGN:
    % SNR_eff = ||h||^2 * Es / N0_eff  = Es / N0  -> N0_eff = ||h||^2 * N0
    sigma_isi = sqrt(h_energy * N0 / 2);

    % Random bits
    N_bits_qpsk = floor(N_bits/2)*2;
    b = randi([0 1], 1, N_bits_qpsk);
    N_sym = N_bits_qpsk/2;

    % Bits -> QPSK symbols
    b_mat = reshape(b, 2, N_sym).';
    x = zeros(1, N_sym);
    for k = 1:N_sym
        x(k) = qpsk_mod(b_mat(k,:), A);
    end

    %% ---------- AWGN BASELINE ----------
    n_awgn = sigma_awgn * (randn(size(x)) + 1j*randn(size(x)));
    y_awgn = x + n_awgn;

    b_hat_mat = zeros(N_sym,2);
    for k = 1:N_sym
        b_hat_mat(k,:) = qpsk_demod(y_awgn(k));
    end
    b_hat = reshape(b_hat_mat.',1,[]);
    BER_AWGN(idx) = mean(b ~= b_hat);

    %% ---------- ISI CHANNEL ----------
    y = zeros(1, N_sym);
    n_isi = sigma_isi * (randn(1,N_sym) + 1j*randn(1,N_sym));

    % y_1 = 2 x_1 + n_1 (no previous symbol)
    y(1) = 2*x(1) + n_isi(1);
    % y_k = 2 x_k + x_{k-1} + n_k
    for k = 2:N_sym
        y(k) = 2*x(k) + x(k-1) + n_isi(k);
    end

    %% ---------- VITERBI MLSE ----------
    numStates  = 4;
    pathMetric = zeros(1, numStates);
    prevState  = zeros(numStates, N_sym);

    % Initialization at k = 1
    % We don't know x_0, so for each candidate x_1=S(s_curr)
    % we pick the best match over all possible x_0 in S.
    for s_curr = 1:numStates
        x_curr = S(s_curr);
        best_metric = inf;

        for s_prev = 1:numStates
            x_prev = S(s_prev);
            y_hat = 2*x_curr + x_prev;
            m = abs(y(1) - y_hat)^2;
            if m < best_metric
                best_metric = m;
            end
        end

        pathMetric(s_curr) = best_metric;
        prevState(s_curr,1) = 0;   % no real previous at k=1
    end

    % Recursion for k >= 2
    for k = 2:N_sym
        path_new = inf(1, numStates);
        prev_new = zeros(1, numStates);

        for s_curr = 1:numStates
            x_curr = S(s_curr);

            best_metric = inf;
            best_prev   = 1;

            for s_prev = 1:numStates
                x_prev = S(s_prev);

                y_hat = 2*x_curr + x_prev;
                branch_metric = abs(y(k) - y_hat)^2;
                total_metric  = pathMetric(s_prev) + branch_metric;

                if total_metric < best_metric
                    best_metric = total_metric;
                    best_prev   = s_prev;
                end
            end

            path_new(s_curr) = best_metric;
            prev_new(s_curr) = best_prev;
        end

        pathMetric = path_new;
        prevState(:,k) = prev_new.';
    end

    % Traceback
    [~, final_state] = min(pathMetric);
    state_seq = zeros(1, N_sym);
    state_seq(N_sym) = final_state;

    for k = N_sym:-1:2
        state_seq(k-1) = prevState(state_seq(k), k);
    end

    x_hat = S(state_seq);

    % Demodulate detected symbols
    b_hat_mat = zeros(N_sym, 2);
    for k = 1:N_sym
        b_hat_mat(k,:) = qpsk_demod(x_hat(k));
    end
    b_hat = reshape(b_hat_mat.', 1, []);

    BER_ISI_VIT(idx) = mean(b ~= b_hat);

    fprintf('Eb/N0 %2d dB | AWGN = %.5f | ISI+Viterbi = %.5f\n', ...
        EbN0_dB(idx), BER_AWGN(idx), BER_ISI_VIT(idx));
end

%% ---------- PLOT ----------
figure('Position', [100 100 900 600]);
semilogy(EbN0_dB, BER_AWGN, 'b-o', 'LineWidth', 2, 'MarkerSize', 7, ...
    'DisplayName', 'AWGN (Simulated)'); hold on;

semilogy(EbN0_dB, BER_ISI_VIT, 'r-s', 'LineWidth', 2, 'MarkerSize', 7, ...
    'DisplayName', 'ISI + Viterbi MLSE');

% Theoretical AWGN curve using Q(x) = 0.5 erfc(x/sqrt(2))
Q = @(x) 0.5 * erfc(x ./ sqrt(2));
BER_theory = Q(sqrt(2*10.^(EbN0_dB/10)));

semilogy(EbN0_dB, BER_theory, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'AWGN (Theory)');

grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 13);
ylabel('Bit Error Rate (BER)', 'FontSize', 13);
title('QPSK over ISI Channel: y_i = 2x_i + x_{i-1} + n_i', 'FontSize', 15);
legend('Location', 'southwest', 'FontSize', 11);
ylim([1e-5 1]);

fprintf('\n=== Simulation Complete ===\n');
