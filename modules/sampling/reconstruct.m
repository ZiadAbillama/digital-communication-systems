function xrcon = reconstruct(t, x_sample, fs)
% RECONSTRUCT - Reconstructs a signal from its samples using sinc interpolation.
%
% Inputs:
%   t         : Continuous-time vector for reconstruction
%   x_sample  : Sampled signal values
%   fs        : Sampling frequency
%
% Output:
%   xrcon     : Reconstructed signal

    Ts = 1 / fs;                           % Sampling period
    n = 0:length(x_sample)-1;              % Sample indices
    t_sample = n * Ts;                     % Sample times

    % Sinc reconstruction
    xrcon = zeros(size(t));
    for k = 1:length(x_sample)
        xrcon = xrcon + x_sample(k) * sinc(fs * (t - t_sample(k)));
    end
end
