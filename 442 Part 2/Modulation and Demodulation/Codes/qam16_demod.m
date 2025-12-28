function b = qam16_demod(x)
% 16-QAM demodulation (Gray decoding)

    % Undo normalization
    I = real(x) * sqrt(10);
    Q = imag(x) * sqrt(10);

    % Decide amplitude
    if I < -2, I_val = -3;
    elseif I < 0, I_val = -1;
    elseif I < 2, I_val = 1;
    else, I_val = 3;
    end

    if Q < -2, Q_val = -3;
    elseif Q < 0, Q_val = -1;
    elseif Q < 2, Q_val = 1;
    else, Q_val = 3;
    end

    % Convert amplitude to natural index 0..3
    levels = [-3 -1 1 3];
    I_nat = find(levels == I_val) - 1; % this line gives back index where == is located
    Q_nat = find(levels == Q_val) - 1;

    % Natural to Gray mapping
    nat2gray = [0 1 3 2];

    I_gray = nat2gray(I_nat + 1);
    Q_gray = nat2gray(Q_nat + 1);

    % Convert Gray code to bits
    I_bits = [bitget(I_gray,2) bitget(I_gray,1)];
    Q_bits = [bitget(Q_gray,2) bitget(Q_gray,1)];

    b = [I_bits Q_bits];
end
