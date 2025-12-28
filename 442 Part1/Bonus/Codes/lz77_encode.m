function [encoded, CR, bps] = lz77_encode(data, window_size, lookahead_size)
% Standard LZ77 encoder (guaranteed to compress repeating patterns)

data = data(:)';          
n = length(data);
encoded = [];

i = 1;
while i <= n
    max_len = 0;
    max_offset = 0;

    win_start = max(1, i - window_size);
    window = data(win_start:i-1);

    % Search longest match in window
    for offset = 1:length(window)
        pos = i - offset;
        len = 0;

        while i+len <= n && data(pos+len) == data(i+len)
            len = len + 1;
            if len == lookahead_size
                break;
            end
        end

        if len > max_len
            max_len = len;
            max_offset = offset;
        end
    end

    % Emit token
    if max_len >= 2
        encoded = [encoded; 1, max_offset, max_len];
        i = i + max_len;
    else
        encoded = [encoded; 0, 0, data(i)];
        i = i + 1;
    end
end

% ----- Bit accounting -----
A = unique(data);
bSym = max(1, ceil(log2(numel(A))));
bOff = max(1, ceil(log2(window_size)));
bLen = max(1, ceil(log2(lookahead_size)));

flags = encoded(:,1);
L = sum(flags==0);
M = sum(flags==1);

bits = ... 
    size(encoded,1) * 1 + ...     % flag
    L * bSym + ...                % literal symbol
    M * (bOff + bLen);            % offset + length

uncompressed = n * bSym;

bps = bits / n;
CR = uncompressed / bits;

end
