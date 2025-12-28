
function [xq, thr, lvl] = lloyd_max(x, M, tol, max_iter)
% Function: lloyd_max.m
% Description:
%   Implements the Lloyd–Max optimal quantizer for minimizing mean squared error.
%
% Inputs:
%   x         - input signal (vector)
%   M         - number of quantization levels
%   tol       - convergence tolerance (optional, default = 1e-6)
%   max_iter  - maximum iterations (optional, default = 100)
%
% Outputs:
%   xq  - quantized signal
%   thr - decision thresholds
%   lvl - quantization levels (centroids)

if nargin < 3, tol = 1e-6; end
if nargin < 4, max_iter = 100; end

x = x(:); % ensure column vector

% Step 1: Initialize thresholds and levels (uniformly)
xmin = min(x);
xmax = max(x);
thr = linspace(xmin, xmax, M+1);
thr = thr(2:end-1);
lvl = linspace(xmin, xmax, M);

prev_lvl = lvl;
err = inf;
iter = 0;

% Step 2: Iterative Lloyd–Max optimization
while err > tol && iter < max_iter
    iter = iter + 1;
    
    % Update thresholds as midpoints between levels
    thr = 0.5 * (lvl(1:end-1) + lvl(2:end));
    
    % Assign each sample to nearest level
    xq = quan(x, thr, lvl);
    
    % Update levels as centroids (mean of assigned samples)
    for k = 1:length(lvl)
        if k == 1
            idx = x <= thr(1);
        elseif k == length(lvl)
            idx = x > thr(end);
        else
            idx = x > thr(k-1) & x <= thr(k);
        end
        if any(idx)
            lvl(k) = mean(x(idx));
        end
    end
    
    % Check convergence
    err = max(abs(lvl - prev_lvl));
    prev_lvl = lvl;
end

% Final quantization
xq = quan(x, thr, lvl);

end
