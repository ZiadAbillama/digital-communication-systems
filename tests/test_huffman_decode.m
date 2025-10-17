% Test: test_huffman_decode.m
% Purpose:
%   Verify Huffman decoding using the dictionary from huffman_encode
%   and ensure lossless reconstruction.

clc; clear; close all;

% Example 1: Simple discrete sequence
fprintf('Example 1: Simple Discrete Source\n');

% Original source
x = [1 1 1 2 2 3 3 3 3 3];

% Build Huffman dictionary
[dict, avglen] = huffman_encode(x);

% Encode using MATLAB's huffmanenco
encoded = huffmanenco(num2cell(x), dict);

% Decode using our function
decoded = huffman_decode(encoded, dict);

% Compute entropy & efficiency
H = entropy(x);
efficiency = (H / avglen) * 100;

% Results
fprintf('Entropy: %.4f bits/symbol\n', H);
fprintf('Average Huffman code length: %.4f bits/symbol\n', avglen);
fprintf('Efficiency: %.2f%%\n', efficiency);
fprintf('Lossless reconstruction: %s\n', string(isequal(x(:), decoded(:))));
disp('Huffman Dictionary:');
disp(dict);



% Example 2: Quantized signal
fprintf('Example 2: Quantized Signal\n');

% Quantize a signal
t = 0:0.001:1;
signal = sin(2*pi*5*t) + 0.3*randn(size(t));
nLevels = 8;
xq = round(rescale(signal, 1, nLevels));

% Build Huffman dictionary for quantized symbols
[dict_q, avglen_q] = huffman_encode(xq);

% Encode using MATLAB's huffmanenco
encoded_q = huffmanenco(num2cell(xq), dict_q);

% Decode using our function
decoded_q = huffman_decode(encoded_q, dict_q);

% Compute entropy & efficiency
Hq = entropy(xq);
efficiency_q = (Hq / avglen_q) * 100;

% Results
fprintf('Entropy: %.4f bits/symbol\n', Hq);
fprintf('Average Huffman code length: %.4f bits/symbol\n', avglen_q);
fprintf('Efficiency: %.2f%%\n', efficiency_q);
fprintf('Lossless reconstruction: %s\n', string(isequal(xq(:), decoded_q(:))));
disp('Huffman Dictionary (Quantized Signal):');
disp(dict_q);

disp('✅ test_huffman_decode.m executed successfully.');
