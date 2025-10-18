clc; clear; close all;
disp('Test: LZ77 Encoding');

x = [1 1 1 2 2 3 3 3 3 3];
[stream, CR, bits] = lz77_encode(x, 8, 4);

disp('[offset, length, nextSymbol]:');
disp(stream);
fprintf('Compression Ratio = %.2f\n', CR);
fprintf('Avg bits/symbol   = %.2f\n', bits);
disp('✅ test_lz77_encode.m executed successfully.');
