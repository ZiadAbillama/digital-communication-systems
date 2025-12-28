function b_hat = repetition_decode(r, L)
% Repetition decoder (majority vote)
% r : received bits after channel (0/1 values or soft values)
% L : repetition factor used in the encoder
% b_hat : decoded bit vector

    % make sure r is a row vector
    r = r(:)';

    % number of symbols
    N = length(r) / L;

    % pre-allocate output bits
    b_hat = zeros(1, N);

    % decode each group of L bits
    for i = 1:N
        
        % extract the block of L received bits
        block = r((i-1)*L + 1 : i*L);

        % majority vote:
        % if more than half are >= 0.5 → bit is 1
        if sum(block >= 0.5) > L/2 % block>= 0.5 returns 1 if T 0 if F for each bit
            b_hat(i) = 1;
        else
            b_hat(i) = 0;
        end
    end

end
