% Test: test_huffman_encode.m
% Purpose:
%   Verify Huffman encoding and compare to entropy.

clc; clear; close all;

% Example 1: simple discrete sequence
fprintf('Example 1: Discrete source sequence');
x = [1 1 1 2 2 3 3 3 3 3];
H = entropy(x);
[dict, avglen] = huffman_encode(x);

fprintf('Entropy (bits per symbol): %.4f\n', H);
fprintf('Average Huffman code length: %.4f\n', avglen);
fprintf('Efficiency: %.2f%%\n', (H/avglen)*100);

% Display dictionary
disp('Huffman Dictionary:');
disp(dict);

% Example 2: quantized signal
fprintf('Example 2: Quantized signal');
t = 0:0.001:1;
signal = sin(2*pi*5*t) + 0.3*randn(size(t));
nLevels = 8;
xq = round(rescale(signal, 1, nLevels));
Hq = entropy(xq);
[dict_q, avglen_q] = huffman_encode(xq);
fprintf('Entropy (bits/symbol): %.4f\n', Hq);
fprintf('Average Huffman code length: %.4f\n', avglen_q);
fprintf('Efficiency: %.2f%%\n', (Hq/avglen_q)*100);

% Display dictionary
disp('Huffman Dictionary (Quantized Signal):');
disp(dict_q);