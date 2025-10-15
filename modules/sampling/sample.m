function [t_sample, x_sample] = sample(t, xt, fs)
% sample - Samples a continuous-time signal at a specified rate
% Inputs:
%   t  - original time vector
%   xt - original signal values at t
%   fs - sampling frequency (Hz)
% Outputs:
%   t_sample - sampled time points
%   x_sample - sampled signal values

    Ts = 1/fs;                     % Sampling period
    t_sample = t(1) : Ts : t(end);   % New sampled time points
    x_sample = interp1(t, xt, t_sample, 'linear');  % Interpolate to find sample values
end
