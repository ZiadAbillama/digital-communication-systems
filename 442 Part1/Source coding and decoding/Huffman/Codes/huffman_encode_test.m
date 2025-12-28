% Example 1 – Discrete source
fprintf('\n Example 1: Simple discrete sequence\n');
x = [1 1 1 2 2 3 3 3 3 3];

% PMF and entropy, also computed in function but here for clarity
% during testing, the other serves as meta data
symbols = unique(x); %finds unique symobls in x
counts = histcounts(categorical(x), categorical(symbols)); %Counts how often each occurs
p = counts / numel(x); %probability
H = -sum(p .* log2(p + eps)); %+eps is safety from 0

% Encode using manual Huffman
[encoded_bits, dict, avglen] = huffman_encode(x);
disp('Encoded bitstream:');
disp(encoded_bits);

% Fixed-length baseline
L_fixed = ceil(log2(dict.M));

fprintf('Entropy (bits/symbol): %.4f\n', H);
fprintf('Average Huffman length: %.4f\n', avglen);
fprintf('Fixed-length baseline: %d bits\n', L_fixed);
fprintf('Efficiency: %.2f%%\n', (H/avglen)*100);

disp('Huffman Dictionary:');
disp(dict);

% Example 2 – Quantized signal from uniform quantizer
fprintf('Example 2: Unifromly Quantized signal\n');
fs = 1000; t = 0:1/fs:1;
signal = sin(2*pi*5*t) + 0.3*randn(size(t));

M = 8;
xmin = min(signal); xmax = max(signal);
thr = linspace(xmin, xmax, M+1); thr = thr(2:end-1);
lvl = linspace(xmin+(xmax-xmin)/(2*M), xmax-(xmax-xmin)/(2*M), M);
xq = quan(signal, thr, lvl);

% Compute entropy manually
symbols_q = unique(xq);
counts_q = histcounts(categorical(xq), categorical(symbols_q));
p_q = counts_q / numel(xq);
Hq = -sum(p_q .* log2(p_q + eps));

% Encode
[encoded_bits_q, dict_q, avglen_q] = huffman_encode(xq);
disp('Encoded Huffman bitstream (first 50 bits):');
disp(encoded_bits_q(1:min(50, length(encoded_bits_q))));
L_fixed_q = ceil(log2(dict_q.M));

fprintf('Entropy (bits/symbol): %.4f\n', Hq);
fprintf('Average Huffman length: %.4f\n', avglen_q);
fprintf('Fixed-length baseline: %d bits\n', L_fixed_q);
fprintf('Efficiency: %.2f%%\n', (Hq/avglen_q)*100);

disp('Huffman Dictionary (Quantized Signal):');
disp(dict_q);

% Example 3 – Lloyd–Max quantizer
fprintf('Example 3: Lloyd–Max quantized signal\n');

% Perform Lloyd–Max quantization using the same signal and M
[xq_lloyd, thr_lloyd, lvl_lloyd] = lloyd_max(signal, M);

% Compute PMF and entropy for Lloyd–Max levels
symbols_l = unique(xq_lloyd);
counts_l = histcounts(categorical(xq_lloyd), categorical(symbols_l));
p_l = counts_l / numel(xq_lloyd);
Hl = -sum(p_l .* log2(p_l + eps));

% Huffman encode using the Lloyd–Max quantized sequence
[encoded_bits_l, dict_l, avglen_l] = huffman_encode(xq_lloyd);

% Display encoded bitstream and metrics
disp('Encoded Huffman bitstream (first 50 bits):');
disp(encoded_bits_l(1:min(50, length(encoded_bits_l))));
L_fixed_l = ceil(log2(dict_l.M));
fprintf('Entropy (bits/symbol): %.4f\n', Hl);
fprintf('Average Huffman length: %.4f\n', avglen_l);
fprintf('Fixed-length baseline: %d bits\n', L_fixed_l);
fprintf('Efficiency: %.2f%%\n', (Hl/avglen_l)*100);

% Display Huffman dictionary for Lloyd–Max quantization
disp('Huffman Dictionary (Lloyd–Max Quantized Signal):');
disp(dict_l);
