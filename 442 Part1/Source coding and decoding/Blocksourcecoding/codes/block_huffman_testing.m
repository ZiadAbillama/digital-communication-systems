% Test: test_block_huffman.m
disp('--- Test: Block Huffman Coding ---');

% Example A: simple discrete source 
x = [1 1 1 2 2 3 3 3 3 3];
for k = [1 2 3 4]
    [~, dictA, avglenA, effA, HkA] = block_huffman(x, k);
    fprintf('\n[Discrete] k=%d | Hk/k=%.4f | Lbar=%.4f | Eff=%.2f%% | |A_k|=%d\n', ...
        k, HkA, avglenA, effA, dictA.M);
end

%  Example B: uniform quantizer 
fs = 1000; t = 0:1/fs:1;
signal = sin(2*pi*5*t) + 0.3*randn(size(t));
M = 8;

xmin = min(signal); xmax = max(signal);
thr = linspace(xmin, xmax, M+1); thr = thr(2:end-1);
lvl = linspace(xmin+(xmax-xmin)/(2*M), xmax-(xmax-xmin)/(2*M), M);
xq_uni = quan(signal, thr, lvl);

for k = [2 3 4]
    [~, dictU, avglenU, effU, HkU] = block_huffman(xq_uni, k);
    fprintf('\n[Uniform] k=%d | Hk/k=%.4f | Lbar=%.4f | Eff=%.2f%% | |A_k|=%d\n', ...
        k, HkU, avglenU, effU, dictU.M);
end

%Example C: Lloyd–Max quantizer 
[xq_lm, thr_lm, lvl_lm] = lloyd_max(signal, M);
for k = [2 3 4]
    [~, dictL, avglenL, effL, HkL] = block_huffman(xq_lm, k);
    fprintf('\n[Lloyd–Max] k=%d | Hk/k=%.4f | Lbar=%.4f | Eff=%.2f%% | |A_k|=%d\n', ...
        k, HkL, avglenL, effL, dictL.M);
end
