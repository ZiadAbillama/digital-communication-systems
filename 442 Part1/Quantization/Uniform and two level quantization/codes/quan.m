function xq = quan(x, thr, lvl)
% Quantize samples x using decision thresholds thr and reconstruction levels lvl.
% It maps the values x to the representation points lvl
% Inputs:
%   x   : row/col vector 
%   thr : Threshold (Boundaries) strictly increasing
%   lvl : Representation points
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
    idx = sum(xv > thr, 2) + 1;   % gives which region the points falls in

    % Map to levels and reshape
    xq = reshape(lvl(idx), sz);
end
