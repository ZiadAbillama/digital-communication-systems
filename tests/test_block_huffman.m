clc; clear; close all;
disp('--- Test: Block Huffman Coding ---');

% Example discrete source
x = [1 1 1 2 2 3 3 3 3 3];

for k = [1 2 3 4]
    [dict, avglen, eff, Hk] = block_huffman(x, k);
    fprintf('\nBlock length k=%d:\n', k);
    fprintf('Entropy rate per symbol Hk/k = %.4f bits/symbol\n', Hk);
    fprintf('Average code length Lbar = %.4f bits/symbol\n', avglen);
    fprintf('Efficiency = %.2f%%\n', eff);
end

disp('✅ test_block_huffman.m executed successfully.');
