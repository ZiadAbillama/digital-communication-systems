function xrcon = reconstruct(t, x_sample, fs)
% RECONSTRUCT - Rebuilds a continuous-time signal from its discrete samples
% Inputs:
%   t         : Continuous-time vector where the signal is reconstructed
%   x_sample  : Discrete sampled values of the original signal
%   fs        : Sampling frequency (in Hz)
%
% Output:
%   xrcon     : Reconstructed signal evaluated at times in t

    % Makes sure inputs are row vectors
    t        = t(:).';
    x_sample = x_sample(:).';

    Ts = 1/fs;                          % Sampling period
    n  = 0:length(x_sample)-1;          % Sample indices 
    t_sample = t(1) + n*Ts;             % Actual sample times aligned with sample() output

    xrcon = zeros(size(t));             % Initialize reconstructed signal with zeros

    for k = 1:length(x_sample)
        % Compute time difference between t and current sample
        x_diff = (t - t_sample(k)) / Ts;

        % Create custom sinc
        sinc_val = ones(size(x_diff));      % Start with 1s (handles x = 0 safely)
        nz = (x_diff ~= 0);                 % Identify non-zero indices
        sinc_val(nz) = sin(pi * x_diff(nz)) ./ (pi * x_diff(nz));

        % Add this sample's sinc contribution to the reconstructed signal
        xrcon = xrcon + x_sample(k) * sinc_val;
    end
end
