function [t_sample, x_sample] = sample(t, xt, fs)
% Function: sample.m
% Description:
%   Samples a continuous-time signal x(t) at a given sampling frequency fs.
%
% Inputs:
%   t  - time vector of the continuous signal
%   xt - signal values corresponding to t
%   fs - sampling frequency (Hz)
%
% Outputs:
%   t_sample - time instances where signal is sampled
%   x_sample - sampled signal values

% Ensure t and xt are row vectors
t = t(:).';
xt = xt(:).';

% Sampling period
Ts = 1/fs;

% Create sampling instants within the range of t
t_sample = t(1):Ts:t(end);

% Interpolate signal values at sampled times
x_sample = interp1(t, xt, t_sample, 'linear');

end
