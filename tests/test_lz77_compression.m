% Test: test_lz77_compression.m
% Purpose:
%   Compare Huffman and LZ77 compression performance
%   on quantized signals with different quantization levels (M).

clc; clear; close all;
disp('--- Test: LZ77 vs Huffman Compression ---');

% Create a representative signal
t = linspace(0, 1, 2000);
x = sin(2*pi*3*t) + 0.3*sin(2*pi*7*t);

% Quantization levels to test
M_values = [4 8 16 32];
CR_huff = zeros(size(M_values));
CR_lz77 = zeros(size(M_values));
bits_huff = zeros(size(M_values));
bits_lz77 = zeros(size(M_values));

for i = 1:length(M_values)
    M = M_values(i);
    % --- Uniform quantization ---
    xmin = min(x); xmax = max(x);
    thr = linspace(xmin, xmax, M+1); thr = thr(2:end-1);
    lvl = linspace(xmin+(xmax-xmin)/(2*M), xmax-(xmax-xmin)/(2*M), M);
    xq = quan(x, thr, lvl);

    % --- Huffman coding ---
    [dict, avglen] = huffman_encode(xq);
    H = entropy(xq);
    Eff_huff = 100 * H / avglen;
    CR_huff(i) = 8 / avglen;         % ratio (assuming 8-bit original)
    bits_huff(i) = avglen;

    % --- LZ77 encoding ---
    [~, CR, bits] = lz77_encode(xq, 32, 16);
    CR_lz77(i) = CR;
    bits_lz77(i) = bits;
end

% --- Print results ---
fprintf('\n%-6s %-12s %-12s %-12s %-12s\n', 'M', 'CR_Huff', 'CR_LZ77', 'Bits_Huff', 'Bits_LZ77');
for i = 1:length(M_values)
    fprintf('%-6d %-12.3f %-12.3f %-12.3f %-12.3f\n', ...
        M_values(i), CR_huff(i), CR_lz77(i), bits_huff(i), bits_lz77(i));
end

% --- Plot compression ratio comparison ---
figure('Name','LZ77_vs_Huffman_Compression','NumberTitle','off');
plot(M_values, CR_huff, 'ro-', 'LineWidth', 1.3); hold on;
plot(M_values, CR_lz77, 'bs-', 'LineWidth', 1.3);
xlabel('Quantization Levels (M)');
ylabel('Compression Ratio (Original bits / Compressed bits)');
title('LZ77 vs Huffman Compression Efficiency');
legend('Huffman','LZ77','Location','best'); grid on;

% --- Save figure automatically ---
resultsDir = fullfile(pwd, 'results', 'figures');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end
saveas(gcf, fullfile(resultsDir, 'LZ77_vs_Huffman_Compression.png'));

disp('✅ test_lz77_compression.m executed successfully.');
