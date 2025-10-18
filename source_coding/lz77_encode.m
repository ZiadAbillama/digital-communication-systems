function [encoded_stream, compression_ratio, avg_bits_per_symbol] = lz77_encode(x, window_size, lookahead_size)
% Function: lz77_encode.m
% Description:
%   Simple LZ77 dictionary encoder for discrete input sequence x.
%   Returns encoded triplets [offset, length, nextSymbol].
%
% Inputs:
%   x              - input sequence (e.g. quantized samples)
%   window_size    - size of search window (e.g. 20–50)
%   lookahead_size - size of lookahead buffer (e.g. 10)
%
% Outputs:
%   encoded_stream        - Nx3 matrix [offset, match_length, next_symbol]
%   compression_ratio     - ratio of original length / encoded length (approx.)
%   avg_bits_per_symbol   - average bits per symbol (approx.)

x = x(:)';  % ensure row vector
n = length(x);
encoded_stream = [];

i = 1;
while i <= n
    search_start = max(1, i - window_size);
    search_buffer = x(search_start:i-1);
    lookahead_buffer = x(i:min(i+lookahead_size-1, n));

    best_offset = 0;
    best_length = 0;

    for offset = 1:length(search_buffer)
        match_length = 0;
        while match_length < length(lookahead_buffer) && ...
              (length(search_buffer) - offset + 1 + match_length) <= length(search_buffer) && ...
              search_buffer(length(search_buffer) - offset + 1 + match_length) == lookahead_buffer(match_length + 1)
            match_length = match_length + 1;
        end

        if match_length > best_length
            best_length = match_length;
            best_offset = offset;
        end
    end

    if i + best_length <= n
        next_symbol = x(i + best_length);
    else
        next_symbol = NaN;
    end

    encoded_stream = [encoded_stream; best_offset, best_length, next_symbol];

    i = i + best_length + 1;
end

% Rough compression metrics
bits_per_token = ceil(log2(window_size)) + ceil(log2(lookahead_size)) + 8; % approximate bits per triplet
total_bits = size(encoded_stream,1) * bits_per_token;
avg_bits_per_symbol = total_bits / n;
compression_ratio = (n * 8) / total_bits;

end