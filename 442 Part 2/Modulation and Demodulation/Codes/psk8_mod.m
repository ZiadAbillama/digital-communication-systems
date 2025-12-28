function x = psk8_mod(b, A)
% 8-PSK modulation
% b: 3-bit input 
% A: amplitude
    b = b(:)';

    % Binary to decimal
    dec = b(1)*4 + b(2)*2 + b(3);

    % 8 equally spaced points starting at 0
    theta = (2*pi/8) * dec;
    x = A * exp(1j * theta);
end