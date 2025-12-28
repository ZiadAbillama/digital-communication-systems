function decoded = huffman_decode(encoded_bits, dict)
% Function: huffman_decode.m
% Inputs:
%   encoded_bits - char row vector of '0'/'1' (from huffman_encode)
%   dict         - struct with fields:
%                  .symbols (1xM), .codes (1xM cell of char)
% Output:
%   decoded      - column vector of reconstructed symbols

    % Ensure a char row vector
    encoded_bits = char(encoded_bits(:).');  

    % Map: code (char) -> symbol
    code2sym = containers.Map(dict.codes, num2cell(dict.symbols));

    decoded = [];           % output sequence
    buffer  = '';           % bit buffer (char)

    for i = 1:numel(encoded_bits)
        buffer(end+1) = encoded_bits(i);   
        if isKey(code2sym, buffer)
            decoded(end+1,1) = code2sym(buffer); 
            buffer = '';                       % reset for next symbol
        end
    end

    if ~isempty(buffer)
        warning('huffman_decode:Incomplete', ...
                'Leftover bits after decoding: possible mismatch of dict/stream.');
    end
end
