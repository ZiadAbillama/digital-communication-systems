function [t_sample, x_sample] = sample(t, xt, fs)
% SAMPLE - Samples a given continuous-time signal at rate fs.
%
% Inputs:
%   t  : Time vector (row vector)
%   xt : Signal samples corresponding to t
%   fs : Sampling frequency (Hz)
%
% Outputs:
%   t_sample : Sampled time indices
%   x_sample : Sampled signal values

    Ts = 1 / fs;                   % Sampling period
    t_min = min(t);
    t_max = max(t);

    % Generate sampling instants
    t_sample = t_min : Ts : t_max;

    % Interpolate x(t) at sampling points
    x_sample = interp1(t, xt, t_sample, 'linear');  % Linear interpolation

end
