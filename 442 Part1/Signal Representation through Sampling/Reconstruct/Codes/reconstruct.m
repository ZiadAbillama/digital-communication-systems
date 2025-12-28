function xrcon = reconstruct(t, x_sample, fs)
% Function: reconstruct.m
% Description:
%   Reconstructs signal from discrete samples
%   using sinc interpolation.
%
% Inputs:
%   t         - time vector for reconstruction (continuous)
%   x_sample  - sampled signal values
%   fs        - sampling frequency (Hz)
%
% Outputs:
%   xrcon     - reconstructed signal corresponding to t

% Ensure column vectors
t = t(:);
x_sample = x_sample(:);

% Sampling period
Ts = 1 / fs;

% Create time indices for samples
n = 0:length(x_sample)-1;
t_sample = n * Ts;

% Preallocate reconstruction output
xrcon = zeros(size(t));

% Ideal sinc reconstruction
for k = 1:length(x_sample)
    % Compute time difference between t and current sample
    x_diff = (t - t_sample(k)) / Ts;

    % Create custom sinc
    sinc_val = ones(size(x_diff));      
    nz = (x_diff ~= 0);                 % Identify non-zero indices
    sinc_val(nz) = sin(pi * x_diff(nz)) ./ (pi * x_diff(nz)); %sin(pi x)/(pi x)

    % Add this sample's sinc contribution to the reconstructed signal
    xrcon = xrcon + x_sample(k) * sinc_val;
end
end