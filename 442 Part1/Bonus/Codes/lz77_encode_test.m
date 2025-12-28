% test_lz77_simple.m
fprintf('=== Simple LZ77 Test ===\n\n');

% Test 1: Simple repetition (should compress well)
x1 = [1 2 3 1 2 3 1 2 3];
fprintf('Test 1: [1 2 3 1 2 3 1 2 3]\n');
[tokens1, CR1, bps1] = lz77_encode(x1, 8, 8);
fprintf('Tokens:\n');
disp(tokens1);
fprintf('Bits/symbol: %.3f, CR: %.3f\n\n', bps1, CR1);

% Test 2: All different (should not compress)
x2 = [1 2 3 4 5 6 7 8 9];
fprintf('Test 2: [1 2 3 4 5 6 7 8 9] (no repetition)\n');
[tokens2, CR2, bps2] = lz77_encode(x2, 8, 8);
fprintf('Tokens:\n');
disp(tokens2);
fprintf('Bits/symbol: %.3f, CR: %.3f\n\n', bps2, CR2);

% Test 3: Long repetition
x3 = repmat([1 2 3], 1, 10);
fprintf('Test 3: [1 2 3] repeated 10 times\n');
[tokens3, CR3, bps3] = lz77_encode(x3, 32, 16);
fprintf('Number of tokens: %d (input length: %d)\n', size(tokens3,1), length(x3));
fprintf('Bits/symbol: %.3f, CR: %.3f\n\n', bps3, CR3);

% Test 4: Self-referential match
x4 = [1 2 1 2 1 2 1 2];
fprintf('Test 4: [1 2 1 2 1 2 1 2] (overlapping pattern)\n');
[tokens4, CR4, bps4] = lz77_encode(x4, 8, 8);
fprintf('Tokens:\n');
disp(tokens4);
fprintf('Bits/symbol: %.3f, CR: %.3f\n', bps4, CR4);
