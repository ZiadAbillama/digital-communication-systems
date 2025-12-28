function xq = quan(x, thr, lvl)
% Quantize samples x using decision thresholds thr and reconstruction levels lvl.
% Works for 2-level and multi-level (uniform or Lloyd–Max) quantizers.
% Inputs:
%   x   : row/col vector (any shape is preserved)
%   thr : (M-1)x1 or 1x(M-1) thresholds, strictly increasing
%   lvl : Mx1 or 1xM levels (numel(lvl) = numel(thr)+1)
% Output:
%   xq  : quantized x (same size as x)

    % keep original shape
    sz  = size(x);
    xv  = x(:);

    thr = thr(:).';         % row vector
    lvl = lvl(:);           % column vector

    % basic checks
    if numel(lvl) ~= numel(thr) + 1
        error('quan:SizeMismatch','numel(lvl) must be numel(thr)+1.');
    end
    if any(diff(thr) <= 0)
        error('quan:ThrOrder','thr must be strictly increasing.');
    end

    % Bin index: count how many thresholds each sample exceeds (ties go to LEFT bin).
    % Regions: (-inf, thr1], (thr1, thr2], ..., (thr_{M-1}, inf)
    idx = sum(xv > thr, 2) + 1;   % gives values in {1,...,M}

    % Map to levels and reshape
    xq = reshape(lvl(idx), sz);
end
