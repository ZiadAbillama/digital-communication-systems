function [xhat, ck] = ffs(xt, t, n, T)
% Function: ffs.m
% Description:
%   Computes the finite Fourier series (FFS) approximation of a
%   time-limited signal x(t).
%
% Inputs:
%   xt - original signal values (row or column vector)
%   t  - time samples corresponding to xt
%   n  - number of Fourier terms on each side (positive/negative)
%   T  - period of the approximation (seconds)
%
% Outputs:
%   xhat - reconstructed/approximated signal from FFS
%   ck   - Fourier coefficients for k = -n:n

% Ensure row vectors
t = t(:).';
xt = xt(:).';

% Initialize coefficient vector
k = -n:n;
ck = zeros(size(k));

% Compute Fourier coefficients numerically using trapezoidal integration
for i = 1:length(k)
    ck(i) = (1/T) * trapz(t, xt .* exp(-1j * 2*pi*k(i) * t / T));
end

% Reconstruct the signal using the computed coefficients
xhat = zeros(size(t));
for i = 1:length(k)
    xhat = xhat + ck(i) * exp(1j * 2*pi*k(i) * t / T);
end

% Keep only the real part (original signal is real-valued)
xhat = real(xhat);

end
