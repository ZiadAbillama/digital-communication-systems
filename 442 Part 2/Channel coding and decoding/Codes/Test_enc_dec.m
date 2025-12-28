% number of bits to test
N = 20000;

% repetition factor
L = 3;

% SNR for the test
EbN0_dB = 2;
EbN0 = 10^(EbN0_dB/10);

A = 1;   % BPSK amplitude

% generate random bits
b = randi([0 1], 1, N);


% Uncoded BPSK

% modulate
x_uncoded = bpsk_mod(b, A);

% Eb for BPSK (1 bit per symbol)
Eb_uncoded = A^2;

% compute noise variance
N0_uncoded = Eb_uncoded / EbN0;
sigma_uncoded = sqrt(N0_uncoded/2);

% add noise
y_uncoded = x_uncoded + sigma_uncoded*(randn(size(x_uncoded)) + 1j*randn(size(x_uncoded)));

% demodulate
b_hat_uncoded = bpsk_demod(real(y_uncoded));

% compute BER
BER_uncoded = mean(b ~= b_hat_uncoded);



%Coded BPSK with repetition L=3

% repeat bits
c = repetition_encode(b, L);

% modulate repeated bits
x_coded = bpsk_mod(c, A);

% Eb for repetition code: same bit energy spread over L reps → Eb same
Eb_coded = A^2;

% noise variance for coded signal (same Eb/N0)
N0_coded = Eb_coded / EbN0;
sigma_coded = sqrt(N0_coded/2);

% add noise
y_coded = x_coded + sigma_coded*(randn(size(x_coded)) + 1j*randn(size(x_coded)));

% demodulate each repeated bit
r = bpsk_demod(real(y_coded));

% decode using majority vote
b_hat_coded = repetition_decode(r, L);

% compute BER
BER_coded = mean(b ~= b_hat_coded);



%Display Results 
fprintf('Uncoded BPSK BER          = %.5f\n', BER_uncoded);
fprintf('Repetition (L=%d) BPSK BER = %.5f\n', L, BER_coded);

if BER_coded == 0
    improvement = Inf;
else
    improvement = BER_uncoded / BER_coded;
end

fprintf('Improvement factor = %.2f× better\n', improvement);
