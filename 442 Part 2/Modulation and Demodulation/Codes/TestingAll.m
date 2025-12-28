% BPSK test
fprintf('Test BPSK\n');
test_bits = [0 1];
for b = test_bits
    x = bpsk_mod(b, 1);      % map bit to +A / -A
    b_hat = bpsk_demod(x);   % detect sign
    fprintf('%d -> %+.1f -> %d\n', b, x, b_hat);
end
fprintf('\n');

% QPSK test 
fprintf('Test QPSK\n');
test_bits_qpsk = [
    0 0;
    0 1;
    1 1;
    1 0
];
for i = 1:size(test_bits_qpsk,1)
    b = test_bits_qpsk(i,:);
    x = qpsk_mod(b, 1);
    b_hat = qpsk_demod(x);
    fprintf('[%d %d] -> (%+.2f,%+.2f) -> [%d %d]\n', ...
        b(1), b(2), real(x), imag(x), b_hat(1), b_hat(2));
end
fprintf('\n');

% 8-PSK test 
fprintf('Test 8-PSK\n');
test_bits_8psk = [
    0 0 0;
    0 0 1;
    0 1 0;
    0 1 1;
    1 0 0;
    1 0 1;
    1 1 0;
    1 1 1
];
for i = 1:size(test_bits_8psk,1)
    b = test_bits_8psk(i,:);
    x = psk8_mod(b, 1);
    b_hat = psk8_demod(x);
    fprintf('[%d %d %d] -> %.1f deg -> [%d %d %d]\n', ...
        b(1), b(2), b(3), angle(x)*180/pi, ...
        b_hat(1), b_hat(2), b_hat(3));
end
fprintf('\n');

% 16-QAM test 
fprintf('Test 16-QAM\n');
test_bits_16qam = [
    0 0 0 0;
    0 0 0 1;
    0 0 1 1;
    0 0 1 0;
    1 1 0 1;
    1 0 1 0;
    1 1 1 1;
    1 0 0 0
];
% Here x is normalized by sqrt(10), so values will be +- 1/(sqrt10) or +-
% 3/(sqrt10), this will be used to fairly compare processes later, but we
% will display below the unnormalized version
for i = 1:size(test_bits_16qam,1)
    b = test_bits_16qam(i,:);
    x = qam16_mod(b, 1);
    b_hat = qam16_demod(x);
    fprintf('[%d %d %d %d] -> (%+.2f,%+.2f) -> [%d %d %d %d]\n', ...
        b(1), b(2), b(3), b(4), real(x), imag(x), ...
        b_hat(1), b_hat(2), b_hat(3), b_hat(4));
 % unnormalized output of x: mult by sqrt(10):
 x_raw = x * sqrt(10);   % undo normalization
 fprintf('    unnormalized: (%+2.0f,%+2.0f)\n', real(x_raw), imag(x_raw));

end

