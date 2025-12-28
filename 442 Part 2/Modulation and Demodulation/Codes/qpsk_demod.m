function b = qpsk_demod(x)
% QPSK demodulation
% x: complex symbol
% b: 2-bit output

re = real(x);
im = imag(x);

if re >= 0 && im >= 0
    b = [1 1];
elseif re < 0 && im >= 0
    b = [1 0];
elseif re < 0 && im < 0
    b = [0 0];
else  % re >= 0 && im < 0
    b = [0 1];
end

end