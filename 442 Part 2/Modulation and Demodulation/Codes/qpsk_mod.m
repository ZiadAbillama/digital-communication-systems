function x = qpsk_mod(b, A)
% QPSK modulation
% b: 2-bit input 
% A: amplitude
% x: complex symbol

b = b(:)';

if b(1) == 1 && b(2) == 1
    x = (A/sqrt(2)) * (1 + 1j);
elseif b(1) == 1 && b(2) == 0
    x = (A/sqrt(2)) * (-1 + 1j);
elseif b(1) == 0 && b(2) == 0
    x = (A/sqrt(2)) * (-1 - 1j);
else  % 01
    x = (A/sqrt(2)) * (1 - 1j);
end

end