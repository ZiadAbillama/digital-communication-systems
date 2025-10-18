function [dict, avglen, eff, Hk_per_symbol] = block_huffman(x, k)
% Function: block_huffman.m
% Description:
%   Performs block-based Huffman coding on sequence x.
%   Groups symbols into blocks of length k, estimates probabilities,
%   builds Huffman codes, and returns the average codeword length per symbol.
%
% Inputs:
%   x - vector of discrete symbols (e.g. quantized levels)
%   k - block length (e.g. 2, 3, or 4)
%
% Outputs:
%   dict            - Huffman dictionary for symbol blocks
%   avglen          - average code length (bits per *symbol*)
%   eff             - efficiency (%) = (H_k/k) / avglen * 100
%   Hk_per_symbol   - entropy rate per symbol (H_k/k)

x = x(:)';   % row vector

% 1. Form symbol blocks of length k
numBlocks = floor(length(x)/k);
if numBlocks < 1
    error('Input too short for block length k=%d', k);
end
blocks = reshape(x(1:numBlocks*k), k, numBlocks)';
% Represent each block as a string key
blockKeys = cell(numBlocks,1);
for i = 1:numBlocks
    blockKeys{i} = mat2str(blocks(i,:));
end

% 2. Compute empirical probabilities
[uniqueBlocks, ~, idx] = unique(blockKeys);
counts = accumarray(idx, 1);
p = counts / sum(counts);

% 3. Build Huffman dictionary
dict = huffmandict(uniqueBlocks, p);

% 4. Average codeword length (bits per symbol)
code_lengths = cellfun(@length, dict(:,2));
Lbar_per_block = sum(p .* code_lengths);   % bits per *block*
avglen = Lbar_per_block / k;               % bits per *symbol*

% 5. Compute block entropy and efficiency
Hk = -sum(p .* log2(p));                   % bits per *block*
Hk_per_symbol = Hk / k;                    % bits per *symbol*
eff = 100 * (Hk_per_symbol / avglen);

end
