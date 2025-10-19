function [xhat, ck] = ffs(xt, t, n, T)
% Finite Fourier Series over one time-limited chunk of length T.
% Inputs:
%   xt - original signal values
%   t  - corresponding time samples (should cover at least one period)
%   n  - number of harmonics on each side
%   T  - period of the approximation
%
% Outputs:
%   xhat - reconstructed (approximate) signal
%   ck   - Fourier coefficients for k = -n:n

% Ensure row vectors
t  = t(:).';
xt = xt(:).';

% Restrict to one analysis window [-T/2, T/2]
mask = (t >= -T/2) & (t <= T/2);
tc   = t(mask);
xc   = xt(mask);

if numel(tc) < 2
    error('ffs:NotEnoughSamples','Need ≥2 samples inside [-T/2, T/2].');
end

k  = -n:n;
ck = zeros(1, numel(k));

% Compute Fourier coefficients numerically
for i = 1:numel(k)
    ck(i) = (1/T) * trapz(tc, xc .* exp(-1j * 2*pi * k(i) * tc / T));
end

% Reconstruct the signal over the full time axis
xhat = zeros(size(t));
for i = 1:numel(k)
    xhat = xhat + ck(i) * exp(1j * 2*pi * k(i) * t / T);
end

xhat = real(xhat); % keep real part

end
