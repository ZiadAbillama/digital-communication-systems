function xq = quan(x, thr, lvl)
% Function: quan.m
% Description:
%   Quantizes input samples based on provided thresholds and levels.
%
% Inputs:
%   x   - input signal (row or column vector)
%   thr - vector of quantization thresholds (length M-1)
%   lvl - vector of quantization levels (length M)
%
% Outputs:
%   xq  - quantized signal (same size as x)

% Ensure column vectors
x  = x(:);
thr = thr(:);
lvl = lvl(:);

% Number of levels
M = length(lvl);

% Initialize output
xq = zeros(size(x));

for i = 1:length(x)
    if x(i) <= thr(1)
        xq(i) = lvl(1);
    elseif x(i) > thr(end)
        xq(i) = lvl(end);
    else
        for k = 1:M-1
            if x(i) > thr(k) && x(i) <= thr(k+1)
                xq(i) = lvl(k+1);
                break;
            end
        end
    end
end

% Reshape to match input
xq = reshape(xq, size(x));

end
