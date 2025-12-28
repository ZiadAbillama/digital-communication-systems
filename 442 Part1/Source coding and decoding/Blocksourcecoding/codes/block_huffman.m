function [encoded_bits, dict, avglen, eff, Hk_per_symbol] = block_huffman(x, k)
% Function: block_huffman.m

    x = x(:).';                          % row vector
    N = numel(x);
    nb = floor(N / k);                   % number of full blocks
    if nb < 1
        error('block_huffman:ShortInput','Sequence too short for k=%d.', k);
    end

    % Blocks as rows (nb x k)
    XB = reshape(x(1:nb*k), k, nb).';

    % Unique block symbols and indices per block
    [U, ~, idx] = unique(XB, 'rows', 'stable');   % U: M x k, idx: nb x 1 in 1..M
    counts = accumarray(idx, 1, [], @sum);
    p = counts / sum(counts);                     % 1 x M (implicitly)

    % Block entropy and per-symbol entropy rate
    Hk = -sum(p .* log2(p + eps));                % bits per block
    Hk_per_symbol = Hk / k;                       % bits per symbol

    % Build Huffman on the block-index stream 
    [bits_block, dict_num, Lbar_block] = huffman_encode(idx); 

    % Per-symbol average length and efficiency
    avglen = Lbar_block / k;                      % bits per symbol
    eff = 100 * (Hk_per_symbol / avglen);

    % Build a block-level dictionary aligned with U (symbols are k-length rows)
    % dict_num.symbols are numeric labels (1..M) in stable order; matches U rows.
    dict.symbols = mat2cell(U, ones(size(U,1),1), size(U,2)); % cell of 1xk row vectors
    dict.p       = p(:).';
    dict.codes   = dict_num.codes;
    dict.len     = dict_num.len;
    dict.entropy = Hk;                             % per block
    dict.M       = numel(p);

    encoded_bits = bits_block;
end
