function [xhat, ck] = ffs(xt, t, n, T)
% ffs - Finite Fourier Series approximation of a time-limited signal
%
% Inputs:
%   xt : row vector of signal values x(t)
%   t  : corresponding time vector
%   n  : number of harmonics on each side (total 2n+1 terms)
%   T  : period used for the Fourier series
%
% Outputs:
%   xhat : reconstructed signal using the finite series
%   ck   : Fourier coefficients c_k (k = -n:n)

    % ensure row vectors
    xt = xt(:).';
    t  = t(:).';
    
    dt = t(2) - t(1);                 % time step for integration
    k = -n:n;                         % harmonic indices
    
    ck = zeros(1, length(k));         % preallocate coefficient array
    
    % compute Fourier coefficients manually
    for idx = 1:length(k)
        ck(idx) = (1/T) * sum(xt .* exp(-1j*2*pi*k(idx)*t/T)) * dt;
    end
    
    % reconstruct signal using the finite sum
    xhat = zeros(size(t));
    for idx = 1:length(k)
        xhat = xhat + ck(idx) * exp(1j*2*pi*k(idx)*t/T);
    end
    
    % real part only (signal is real)
    xhat = real(xhat);
end
