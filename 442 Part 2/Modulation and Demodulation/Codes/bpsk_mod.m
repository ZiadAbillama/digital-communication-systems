function x = bpsk_mod(b, A)
% BPSK modulation: maps bits to +A or -A
%
% Inputs:
%   b : vector of bits (0 or 1)
%   A : amplitude
%
% Output:
%   x : BPSK-modulated symbols

    % Map 0 → -A and 1 → +A
    x = A * (2*b - 1);

end
