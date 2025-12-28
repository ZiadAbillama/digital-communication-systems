function b = psk8_demod(x)
% 8-PSK demodulation

    % Angle in [0, 2*pi]
    theta = angle(x);
    if theta < 0
        theta = theta + 2*pi;
    end

    % Nearest point among 8 constellation positions
    dec = round(theta / (2*pi/8));
    dec = mod(dec, 8);

    % Convert decimal back to 3 bits
    b1 = floor(dec / 4);         % MSB
    b2 = mod(floor(dec / 2), 2); % middle
    b3 = mod(dec, 2);            % LSB

    b = [b1 b2 b3];
end
