% AWGN Channel BER Simulation for BPSK, QPSK, 8-PSK, and 16-QAM
% Tests both coded and uncoded performance
%Simulation Parameters
N_bits = 1e6;           % Number of bits to simulate per SNR
L = 3;                  % Repetition factor for channel coding
SNR_dB = 0:1:10;       % SNR range in dB (Es/N0)
N_SNR = length(SNR_dB);
A = 1;                  % Fixed amplitude (Es = A^2 = 1)

% Pre-allocate BER arrays
BER_BPSK_uncoded  = zeros(1, N_SNR);
BER_BPSK_coded    = zeros(1, N_SNR);
BER_QPSK_uncoded  = zeros(1, N_SNR);
BER_QPSK_coded    = zeros(1, N_SNR);

% Bonus modulations
BER_8PSK_uncoded  = zeros(1, N_SNR);
BER_16QAM_uncoded = zeros(1, N_SNR);

fprintf('Starting AWGN BER Simulation using SNR = 10log10(Es/N0)\n');
fprintf('Simulating %d bits per SNR point\n\n', N_bits);

%BPSK
fprintf('=== BPSK Simulation ===\n');

for idx = 1:N_SNR
    snr_db  = SNR_dB(idx);
    SNR_lin = 10^(snr_db/10);    % Es/N0 linear
    Es = A^2;                    % Symbol energy
    N0 = Es / SNR_lin;           % N0 value
    sigma = sqrt(N0/2);          % Noise variance

    % Generate random bits
    b = randi([0 1], 1, N_bits);

    %Uncoded BPSK 
    x = bpsk_mod(b, A);
    n = sigma * (randn(size(x)) + 1j*randn(size(x)));
    y = x + n;
    b_hat = bpsk_demod(real(y));

    BER_BPSK_uncoded(idx) = mean(b ~= b_hat);

    %Coded BPSK
    c = repetition_encode(b, L);
    x_coded = bpsk_mod(c, A);

    n_coded = sigma * (randn(size(x_coded)) + 1j*randn(size(x_coded)));
    y_coded = x_coded + n_coded;

    c_hat = bpsk_demod(real(y_coded));
    b_hat_decoded = repetition_decode(c_hat, L);

    BER_BPSK_coded(idx) = mean(b ~= b_hat_decoded);

    fprintf('BPSK SNR %2d dB | uncoded = %.6f | coded = %.6f\n', ...
        snr_db, BER_BPSK_uncoded(idx), BER_BPSK_coded(idx));
end


%QPSK 
fprintf('\n=== QPSK Simulation ===\n');

for idx = 1:N_SNR
    snr_db  = SNR_dB(idx);
    SNR_lin = 10^(snr_db/10);
    Es = A^2;
    N0 = Es / SNR_lin;
    sigma = sqrt(N0/2);

    % Generate even number of bits
    N_bits_qpsk = floor(N_bits/2)*2;
    b = randi([0 1], 1, N_bits_qpsk);

    % Group into QPSK bit-pairs
    N_symbols = N_bits_qpsk/2;
    b_matrix = reshape(b, 2, N_symbols).';

    % Uncoded QPSK 
    x = zeros(1, N_symbols);
    for i = 1:N_symbols
        x(i) = qpsk_mod(b_matrix(i,:), A);
    end

    n = sigma * (randn(size(x)) + 1j*randn(size(x)));
    y = x + n;

    c_hat = zeros(1, N_bits_qpsk);
    for i = 1:N_symbols
        pair_hat = qpsk_demod(y(i));
        c_hat(2*i-1:2*i) = pair_hat;
    end

    BER_QPSK_uncoded(idx) = mean(b ~= c_hat);

    % Coded QPSK (Symbol-level repetition)
    % Repeat each QPSK symbol L times in SEQUENCE
    b_coded_matrix = kron(b_matrix, ones(L,1)); 

    N_symbols_coded = size(b_coded_matrix, 1);

    % Modulate repeated symbols
    x_coded = zeros(1, N_symbols_coded);
    for i = 1:N_symbols_coded
        x_coded(i) = qpsk_mod(b_coded_matrix(i,:), A);
    end

    % AWGN
    n_coded = sigma * (randn(size(x_coded)) + 1j*randn(size(x_coded)));
    y_coded = x_coded + n_coded;

    % Demod
    c_hat_matrix = zeros(N_symbols_coded, 2);
    for i = 1:N_symbols_coded
        c_hat_matrix(i,:) = qpsk_demod(y_coded(i));
    end

    % Reshape into (L repetitions) × (N_symbols) × (2 bits)
    c_hat_reshaped = reshape(c_hat_matrix, L, N_symbols, 2);

    % Majority vote per symbol
    b_hat_matrix = zeros(N_symbols, 2);
    for s = 1:N_symbols
        b_hat_matrix(s,1) = sum(c_hat_reshaped(:,s,1)) > L/2;
        b_hat_matrix(s,2) = sum(c_hat_reshaped(:,s,2)) > L/2;
    end

    % Back to bitstream
    b_hat = reshape(b_hat_matrix.', 1, []);

    BER_QPSK_coded(idx) = mean(b ~= b_hat);

    fprintf('QPSK SNR %2d dB | uncoded = %.6f | coded = %.6f\n', ...
        snr_db, BER_QPSK_uncoded(idx), BER_QPSK_coded(idx));
end


%Extra: 8-PSK 
fprintf('\n=== Extra: 8-PSK ===\n');

for idx = 1:N_SNR
    snr_db  = SNR_dB(idx);
    SNR_lin = 10^(snr_db/10);
    Es = A^2;
    N0 = Es / SNR_lin;
    sigma = sqrt(N0/2);

    % Bits must be multiple of 3
    N_bits_8 = floor(N_bits/3)*3;
    b = randi([0 1], 1, N_bits_8);

    N_symbols = N_bits_8/3;
    b3 = reshape(b, 3, N_symbols).';

    x = zeros(1, N_symbols);
    for i = 1:N_symbols
        x(i) = psk8_mod(b3(i,:), A);
    end

    n = sigma*(randn(size(x)) + 1j*randn(size(x)));
    y = x + n;

    out = zeros(size(b3));
    for i = 1:N_symbols
        out(i,:) = psk8_demod(y(i));
    end

    BER_8PSK_uncoded(idx) = mean(b3(:) ~= out(:));
end


% Extra: 16-QAM
fprintf('\n=== Extra: 16-QAM ===\n');

for idx = 1:N_SNR
    snr_db  = SNR_dB(idx);
    SNR_lin = 10^(snr_db/10);
    Es = A^2;
    N0 = Es / SNR_lin;
    sigma = sqrt(N0/2);

    % Bits must be multiple of 4
    N_bits_16 = floor(N_bits/4)*4;
    b = randi([0 1], 1, N_bits_16);

    N_symbols = N_bits_16/4;
    b4 = reshape(b, 4, N_symbols).';

    x = zeros(1, N_symbols);
    for i = 1:N_symbols
        x(i) = qam16_mod(b4(i,:), A);
    end

    n = sigma*(randn(size(x)) + 1j*randn(size(x)));
    y = x + n;

    out = zeros(size(b4));
    for i = 1:N_symbols
        out(i,:) = qam16_demod(y(i));
    end

    BER_16QAM_uncoded(idx) = mean(b4(:) ~= out(:));
end


%PLOTTING 

figure('Position', [100 100 1200 500]);

subplot(1,2,1);
semilogy(SNR_dB, BER_BPSK_uncoded, 'b-o','LineWidth',1.5); hold on;
semilogy(SNR_dB, BER_BPSK_coded,   'r-s','LineWidth',1.5);
grid on;
xlabel('SNR = 10log_{10}(E_s/N_0) [dB]');
ylabel('BER');
title('BPSK Performance');
legend('Uncoded','Coded (L=3)');
ylim([1e-6 1]);

subplot(1,2,2);
semilogy(SNR_dB, BER_QPSK_uncoded, 'b-o','LineWidth',1.5); hold on;
semilogy(SNR_dB, BER_QPSK_coded,   'r-s','LineWidth',1.5);
grid on;
xlabel('SNR = 10log_{10}(E_s/N_0) [dB]');
ylabel('BER');
title('QPSK Performance');
legend('Uncoded','Coded (L=3)');
ylim([1e-6 1]);

figure;
semilogy(SNR_dB, BER_8PSK_uncoded,'g-o','LineWidth',1.5); hold on;
semilogy(SNR_dB, BER_16QAM_uncoded,'m-s','LineWidth',1.5);
grid on;
xlabel('SNR = 10log_{10}(E_s/N_0) [dB]');
ylabel('BER');
title('Extra: 8-PSK & 16-QAM');
legend('8-PSK','16-QAM');
ylim([1e-6 1]);

fprintf('\n=== Simulation Complete ===\n');
