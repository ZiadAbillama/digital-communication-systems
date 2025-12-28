function x = qam16_mod(b, A)
% 16-QAM modulation with Gray coding
% b = [b1 b2 b3 b4]

    b = b(:)';

    % Map bits → Gray PAM-4 index (0..3)
    I_gray = b(1)*2 + b(2);
    Q_gray = b(3)*2 + b(4);

    % Gray to natural mapping
    gray2nat = [0 1 3 2];   % Gray: 00 01 11 10
    I_nat = gray2nat(I_gray + 1);
    Q_nat = gray2nat(Q_gray + 1);

    % Natural to amplitude: {0,1,2,3} → {-3,-1,1,3}
    levels = [-3 -1 1 3];
    I = levels(I_nat + 1);
    Q = levels(Q_nat + 1);

    % Normalize to unit average power because 16QAM has avg pwr=10
    x = (A/sqrt(10)) * (I + 1j*Q);
end
