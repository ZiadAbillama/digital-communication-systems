function [t_rec, x_rec] = reconstruct(t_sample, x_sample, t_rec)
% reconstruct - Reconstruct signal from samples using manual sinc interpolation
% Inputs:
%   t_sample - sampled time points
%   x_sample - sampled signal values
%   t_rec    - desired fine time vector for reconstruction
% Outputs:
%   t_rec    - reconstruction time vector
%   x_rec    - reconstructed signal values

    Ts = t_sample(2) - t_sample(1);  % Sampling interval
    fs = 1 / Ts;                    % Sampling frequency
    x_rec = zeros(size(t_rec));      % Initialize reconstructed signal

    % Define custom sinc function
    sinc_custom = @(x) sin(pi*x)./(pi*x);
    
    for n = 1:length(t_sample)
        x_rec = x_rec + x_sample(n) * sinc_custom(fs * (t_rec - t_sample(n)));
    end
end
