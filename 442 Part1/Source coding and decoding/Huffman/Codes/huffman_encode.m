function [encoded_bits, dict, avglen] = huffman_encode(x)
% Inputs:
%   x - input discrete source sequence (e.g., quantized samples from quan/lloyd_max)
%
% Outputs:
%   encoded_bits - char array of '0'/'1' bits representing the Huffman-encoded stream
%   dict         - structure with fields:
%       .symbols  : 1xM vector of unique symbols
%       .p        : 1xM vector of symbol probabilities
%       .codes    : 1xM cell array of bit strings
%       .len      : 1xM vector of code lengths
%       .entropy  : scalar entropy H(A)
%       .M        : number of unique symbols
%   avglen       - average codeword length (bits/symbol)

    % Preprocessing & PMF
    x = x(:);                     % make column vector
    [symbols,~,idx] = unique(x);  % find unique symbols & their indices
    counts = accumarray(idx,1);   % count each symbol’s occurrences
    p = counts'/numel(x);         % normalize → probabilities
    M = numel(symbols);           % alphabet size


    % Edge case: single symbol 
    if M == 1
        dict.symbols = symbols';
        dict.p = 1;
        dict.codes = {'0'};
        dict.len = 1;
        dict.entropy = 0;
        dict.M = 1;
        encoded_bits = repmat('0', 1, numel(x));
        avglen = 1;
        return;
    end

    % Build Huffman Dictionary 
    dict = local_huffman_dict(symbols', p);

    %  Encode Stream 
    encoded_bits = local_huffman_encode_idx(idx, dict.codes);

    % Average length & entropy
    dict.entropy = -sum(p .* log2(p + eps));
    avglen = numel(encoded_bits) / numel(x);
end


%Helper Functions

function dict = local_huffman_dict(symbols, p)
% Builds a Huffman dictionary from probabilities 
    M = numel(p);
    nodes = num2cell(1:M);   % leaf node indices
    weights = num2cell(p);   % node weights

    % Combine nodes iteratively (Huffman tree)
    while numel(weights) > 1
        [~, order] = sort(cell2mat(weights));
        i = order(1); j = order(2);

        newNode = {nodes{i}, nodes{j}};
        newW = weights{i} + weights{j};

        keep = true(1, numel(weights)); 
        keep([i j]) = false; %set as false to remove them
        nodes = [nodes(keep), {newNode}]; %this does not include the 2 smallest nodes we just got done with
        weights = [weights(keep), {newW}];
    end

    % Depth-first traversal to assign codes
    codes = repmat({''}, 1, M);
    function walk(node, prefix)
        if iscell(node)
            walk(node{1}, [prefix '0']);
            walk(node{2}, [prefix '1']);
        else
            codes{node} = prefix;
        end
    end
    walk(nodes{1}, '');

    % Safety for degenerate cases
    for k = 1:M
        if isempty(codes{k}), codes{k} = '0'; end
    end

    dict.symbols = symbols(:).';
    dict.p = p(:).';
    dict.codes = codes;
    dict.len = cellfun(@length, codes);
    dict.M = M;
end


function bits = local_huffman_encode_idx(idx, codes)
% Encode integer indices (1..M) into a single bitstring '0'/'1'
    out = strings(numel(idx),1);
    for n = 1:numel(idx)
        out(n) = codes{idx(n)};
    end
    bits = char(join(out,""));
end
