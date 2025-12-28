% Test: test_huffman_decode.m
%   Runs three cases:
%     (1) simple discrete source, (2) quantized signal using quan.m, (3) Lloyd–Max quantized signal

% Example 1: Simple discrete sequence 
fprintf('\n--- Example 1: Simple discrete sequence ---\n');

x = [1 1 1 2 2 3 3 3 3 3];

% Encode -> Decode
[encoded_bits, dict, avglen] = huffman_encode(x);
x_hat = huffman_decode(encoded_bits, dict);

% Checks
is_lossless = isequal(x_hat, x(:));
fprintf('Lossless reconstruction: %s\n', string(is_lossless));

fprintf('Encoded bitstream : ');
disp(encoded_bits);

% Baseline vs Huffman (just for context)
L_fixed = ceil(log2(dict.M));
fprintf('Avg Huffman length: %.4f bits/sym | Fixed-length: %d bits/sym | |A|=%d\n', ...
        avglen, L_fixed, dict.M);

assert(is_lossless, 'Decoder failed on Example 1 (discrete sequence).');


% Example 2: Quantized signal (uniform quantizer)
fprintf('\n--- Example 2: Quantized signal (uniform quantizer) ---\n');

% Generate a test signal
fs = 1000; t = 0:1/fs:1;
signal = sin(2*pi*5*t) + 0.3*randn(size(t));

% Uniform quantizer parameters
M = 8;
xmin = min(signal); xmax = max(signal);
thr = linspace(xmin, xmax, M+1); thr = thr(2:end-1);
lvl = linspace(xmin+(xmax-xmin)/(2*M), xmax-(xmax-xmin)/(2*M), M);

% Quantize with Part 2 function
xq = quan(signal, thr, lvl);

% Encode -> Decode
[encoded_bits_q, dict_q, avglen_q] = huffman_encode(xq);
xq_hat = huffman_decode(encoded_bits_q, dict_q);

% Checks
is_lossless_q = isequal(xq_hat, xq(:));
fprintf('Lossless reconstruction: %s\n', string(is_lossless_q));

% Show a small snippet of the bitstream
fprintf('Encoded bitstream (first 50 bits): ');
disp(encoded_bits_q(1:min(50, length(encoded_bits_q))));

% Baseline vs Huffman and compression ratio (context)
L_fixed_q = ceil(log2(dict_q.M));
total_bits_huff = length(encoded_bits_q);
total_bits_fixed = numel(xq) * L_fixed_q;
compression_gain = 100*(1 - total_bits_huff/total_bits_fixed);

fprintf(['Avg Huffman length: %.4f bits/sym | Fixed-length: %d bits/sym | |A|=%d\n' ...
         'Total bits: Huffman=%d, Fixed=%d | Compression gain vs fixed: %.2f%%\n'], ...
         avglen_q, L_fixed_q, dict_q.M, total_bits_huff, total_bits_fixed, compression_gain);

assert(is_lossless_q, 'Decoder failed on Example 2 (quantized signal).');


% Example 3: Lloyd–Max quantized signal
fprintf('\n--- Example 3: Lloyd–Max quantized signal ---\n');

% Lloyd–Max quantization using same signal and M
[xq_lloyd, thr_lloyd, lvl_lloyd] = lloyd_max(signal, M);

% Encode -> Decode
[encoded_bits_l, dict_l, avglen_l] = huffman_encode(xq_lloyd);
xq_lloyd_hat = huffman_decode(encoded_bits_l, dict_l);

% Check for lossless reconstruction
is_lossless_l = isequal(xq_lloyd_hat, xq_lloyd(:));
fprintf('Lossless reconstruction: %s\n', string(is_lossless_l));

% Show bitstream snippet
fprintf('Encoded bitstream (first 50 bits): ');
disp(encoded_bits_l(1:min(50, length(encoded_bits_l))));

% Compute compression metrics
L_fixed_l = ceil(log2(dict_l.M));
total_bits_huff_l = length(encoded_bits_l);
total_bits_fixed_l = numel(xq_lloyd) * L_fixed_l;
compression_gain_l = 100*(1 - total_bits_huff_l/total_bits_fixed_l);

fprintf(['Avg Huffman length: %.4f bits/sym | Fixed-length: %d bits/sym | |A|=%d\n' ...
         'Total bits: Huffman=%d, Fixed=%d | Compression gain vs fixed: %.2f%%\n'], ...
         avglen_l, L_fixed_l, dict_l.M, total_bits_huff_l, total_bits_fixed_l, compression_gain_l);

assert(is_lossless_l, 'Decoder failed on Example 3 (Lloyd–Max signal).');
