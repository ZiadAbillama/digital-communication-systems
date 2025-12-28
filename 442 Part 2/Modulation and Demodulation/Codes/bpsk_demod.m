function b_hat = bpsk_demod(x)
% BPSK demodulation: decides bits based on sign of x
%
% Input:
%   x : received BPSK symbols
%
% Output:
%   b_hat : detected bits

    %   x >= 0 → 1
    %   x < 0  → 0
    b_hat = x >= 0;

end
