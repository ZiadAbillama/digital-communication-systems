function decoded = huffman_decode(encoded, dict)
% Function: huffman_decode.m
% Description:
%   Decodes a Huffman-encoded bitstream using the given dictionary.
%
% Inputs:
%   encoded_bits - cell array of bit sequences (output from huffman_encode)
%   dict         - Huffman dictionary (symbol, code)
%
% Output:
%   decoded - reconstructed symbol sequence

% Decode using MATLAB's built-in huffmandeco
decoded = huffmandeco(encoded, dict);


end
