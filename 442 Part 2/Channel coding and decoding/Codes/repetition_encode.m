function c = repetition_encode(b, L)
% Repetition encoder
% b : input bit vector (0/1)
% L : repetition factor (ex: 3 means each bit is sent 3 times)
% c : encoded bit vector

    % make sure b is row vector
    b = b(:)';             

    % repeat each bit L times
    % for example b = [1 0], L = 3 → c = [1 1 1 0 0 0]
    c = repelem(b, L);

end
