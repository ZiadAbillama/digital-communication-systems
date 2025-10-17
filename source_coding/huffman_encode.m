function [dict, avglen] = huffman_encode(x)
% Function: huffman_encode.m
% Description:
%   Performs Huffman encoding for a discrete source sequence x.
%
% Inputs:
%   x - input discrete source (vector), e.g., quantized samples
%
% Outputs:
%   encoded_bits - cell array of bit sequences for each symbol
%   dict         - Huffman dictionary (symbol, code)
%   avglen       - average codeword length (bits/symbol)

x = x(:); % ensure column vector

% Step 1: find unique symbols and probabilities
symbols = unique(x);
p = zeros(length(symbols),1);
for i = 1:length(symbols)
    p(i) = sum(x == symbols(i)) / length(x);
end

% Step 2: build Huffman dictionary
dict = huffmandict(num2cell(symbols), p);

% Step 3: compute average codeword length
code_lengths = cellfun(@length, dict(:,2));
avglen = sum(p .* code_lengths);

end
