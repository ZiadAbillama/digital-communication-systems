% QPSK over ISI channel with Viterbi MLSE and Zero-Forcing Equalizer
% Channel model: y_i = h0*x_i + h1*x_{i-1} + n_i
% Simulation Parameters 
N_bits   = 2e5;
EbN0_dB  = 0:2:14;
N_SNR    = length(EbN0_dB);
A        = 1;
Eb       = 0.5;      % 2 bits per symbol
Es       = 1;        % unit-energy constellation

% Channel taps (normalized so |h0|^2 + |h1|^2 = 1)
h_raw = [2 1];
h     = h_raw / norm(h_raw);   % h(1) = h0, h(2) = h1
Lh    = length(h);

%Initialization
BER_AWGN = zeros(1, N_SNR);
BER_VIT  = zeros(1, N_SNR);
BER_ZF   = zeros(1, N_SNR);

% QPSK constellation (Gray)
bits_list = [1 1; 1 0; 0 0; 0 1];
S = zeros(1,4);
for i = 1:4
    S(i) = qpsk_mod(bits_list(i,:), A);
end

% SNR Loop 
for idx = 1:N_SNR
    EbN0_lin = 10^(EbN0_dB(idx) / 10);
    EsN0_lin = 2 * EbN0_lin;          % QPSK: Es = 2Eb
    N0       = Es / EsN0_lin;
    sigma    = sqrt(N0/2);

    % Bits and QPSK symbols
    N_bits_qpsk = floor(N_bits/2)*2;
    b  = randi([0 1], 1, N_bits_qpsk);
    N_sym = N_bits_qpsk / 2;

    b_mat = reshape(b, 2, N_sym).';
    x = zeros(1, N_sym);
    for k = 1:N_sym
        x(k) = qpsk_mod(b_mat(k,:), A);
    end

    %AWGN baseline
    n = sigma*(randn(size(x)) + 1j*randn(size(x)));
    y_awgn = x + n;

    b_hat_mat = zeros(N_sym,2);
    for k = 1:N_sym
        b_hat_mat(k,:) = qpsk_demod(y_awgn(k));
    end
    BER_AWGN(idx) = mean(b ~= reshape(b_hat_mat.',1,[]));

    % ISI channel
    % y_i = h0*x_i + h1*x_{i-1} + n_i, with x_0 = 0
    y     = zeros(1, N_sym);
    n_isi = sigma*(randn(1,N_sym) + 1j*randn(1,N_sym));

    y(1) = h(1)*x(1) + n_isi(1);  % x_0 = 0
    for k = 2:N_sym
        y(k) = h(1)*x(k) + h(2)*x(k-1) + n_isi(k);
    end

    % Viterbi MLSE
    numStates  = 4;
    pathMetric = zeros(1, numStates);
    prevState  = zeros(numStates, N_sym);

    % k = 1: previous symbol x_0 = 0 (known)
    x_prev_known = 0;
    for s_curr = 1:numStates
        x_curr = S(s_curr);
        y_hat  = h(1)*x_curr + h(2)*x_prev_known;
        pathMetric(s_curr) = abs(y(1) - y_hat)^2;
    end

    % k = 2..N_sym: trellis recursion
    for k = 2:N_sym
        path_new = inf(1, numStates);
        prev_new = zeros(1, numStates);

        for s_curr = 1:numStates
            x_curr = S(s_curr);
            best_metric = inf;
            best_prev   = 1;

            for s_prev = 1:numStates
                x_prev = S(s_prev);
                y_hat  = h(1)*x_curr + h(2)*x_prev;

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

        pathMetric    = path_new;
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

    % Demod after MLSE
    b_hat_mat = zeros(N_sym, 2);
    for k = 1:N_sym
        b_hat_mat(k,:) = qpsk_demod(x_hat(k));
    end
    BER_VIT(idx) = mean(b ~= reshape(b_hat_mat.',1,[]));

    % Zero-Forcing equalizer
    Lg    = 5;                % equalizer length
    Lconv = Lg + Lh - 1;

    Hmat = zeros(Lconv, Lg);
    for col = 1:Lg
        Hmat(col:col+Lh-1, col) = h(:);
    end

    d   = ceil(Lconv/2);
    e_d = zeros(Lconv,1);
    e_d(d) = 1;

    gZF = (Hmat' * Hmat) \ (Hmat' * e_d);

    z = conv(y, gZF);
    z = z(d : d+N_sym-1);

    b_hat_mat = zeros(N_sym, 2);
    for k = 1:N_sym
        [~, idx_min] = min(abs(z(k) - S));
        b_hat_mat(k,:) = bits_list(idx_min,:);
    end
    BER_ZF(idx) = mean(b ~= reshape(b_hat_mat.',1,[]));
end

% Plot 
Q = @(x) 0.5 * erfc(x / sqrt(2));
BER_theory = Q(sqrt(2*10.^(EbN0_dB/10)));

figure('Position', [100 100 950 600]);
semilogy(EbN0_dB, BER_AWGN, 'b-o', 'LineWidth', 2, 'MarkerSize', 7); hold on;
semilogy(EbN0_dB, BER_VIT,  'r-s', 'LineWidth', 2, 'MarkerSize', 7);
semilogy(EbN0_dB, BER_ZF,   'm-d', 'LineWidth', 2, 'MarkerSize', 7);
semilogy(EbN0_dB, BER_theory, 'k--', 'LineWidth', 1.5);

xlabel('E_b/N_0 (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('QPSK over ISI Channel with MLSE and ZF Equalization', 'FontSize', 14);
legend('AWGN (Simulated)', 'Viterbi MLSE', 'Zero-Forcing', 'AWGN (Theory)', ...
       'Location', 'southwest');
grid on;
ylim([1e-5 1]);
