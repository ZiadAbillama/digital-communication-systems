fs = 1000; t = 0:1/fs:1;
signal = sin(2*pi*5*t) + 0.3*randn(size(t));
M = 8;

xmin = min(signal); xmax = max(signal);
thr = linspace(xmin, xmax, M+1); thr = thr(2:end-1);
lvl = linspace(xmin+(xmax-xmin)/(2*M), xmax-(xmax-xmin)/(2*M), M);

% Uniform Quantizer Full Chain
xq_uni = quan(signal, thr, lvl);
[bits_uni, dict_uni, Lbar_uni] = huffman_encode(xq_uni);
xq_uni_hat = huffman_decode(bits_uni, dict_uni);         % column

signal_uni_hat = xq_uni_hat.';                           % row

lossless_uni   = isequal(xq_uni_hat, xq_uni(:));
total_bits_uni = length(bits_uni);
throughput_uni = Lbar_uni;                               % bits/symbol
latency_uni    = numel(xq_uni);                          % buffered to learn p
snr_uni        = 10*log10(mean(signal.^2) / mean((signal - signal_uni_hat).^2));

% Lloyd–Max Quantizer Full Chain
[xq_lm, thr_lm, lvl_lm] = lloyd_max(signal, M);
[bits_lm, dict_lm, Lbar_lm] = huffman_encode(xq_lm);
xq_lm_hat   = huffman_decode(bits_lm, dict_lm);
signal_lm_hat = xq_lm_hat.';                             % row

lossless_lm   = isequal(xq_lm_hat, xq_lm(:));
total_bits_lm = length(bits_lm);
throughput_lm = Lbar_lm;
latency_lm    = numel(xq_lm);
snr_lm        = 10*log10(mean(signal.^2) / mean((signal - signal_lm_hat).^2));

% Block Coding (for comparison)
k = 2;
[~, dictU2, LbarU2, effU2, HkU2] = block_huffman(xq_uni, k);
[~, dictL2, LbarL2, effL2, HkL2] = block_huffman(xq_lm,  k);

% Report 
fprintf('\n=== Integration (Uncoded Channel) ===\n');

fprintf('\n[Uniform] symbol-wise Huffman:\n');
fprintf('Lossless: %d | Total bits: %d | bits/symbol: %.4f | latency(syms): %d\n', ...
    lossless_uni, total_bits_uni, throughput_uni, latency_uni);
fprintf('SNR: %.2f dB | |A|=%d\n', snr_uni, dict_uni.M);

fprintf('\n[Lloyd–Max] symbol-wise Huffman:\n');
fprintf('Lossless: %d | Total bits: %d | bits/symbol: %.4f | latency(syms): %d\n', ...
    lossless_lm, total_bits_lm, throughput_lm, latency_lm);
fprintf('SNR: %.2f dB | |A|=%d\n', snr_lm, dict_lm.M);

fprintf('\n[Comparison] (same seed/noise, same M, block length = 1):\n');
fprintf('Uniform:  bits/sym=%.4f, Total bits=%d\n', throughput_uni, total_bits_uni);
fprintf('Lloyd–Max: bits/sym=%.4f, Total bits=%d\n', throughput_lm, total_bits_lm);
fprintf('Rate ratio (LM/Uniform): %.3f\n', throughput_lm/throughput_uni);

fprintf('\n[Block Huffman (k=%d)]:\n', k);
fprintf('Uniform:  Hk/k=%.4f | bits/sym=%.4f | Eff=%.2f%% | |A_k|=%d\n', ...
    HkU2, LbarU2, effU2, dictU2.M);
fprintf('Lloyd–Max: Hk/k=%.4f | bits/sym=%.4f | Eff=%.2f%% | |A_k|=%d\n', ...
    HkL2, LbarL2, effL2, dictL2.M);
