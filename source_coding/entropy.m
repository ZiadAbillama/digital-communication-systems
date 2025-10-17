function H = entropy(x)
% Function: entropy.m
% Description:
%   Computes the entropy of a discrete source sequence x.
%
% Input:
%   x - Discrete signal (vector), e.g., quantized output levels
%
% Output:
%   H - Entropy in bits per symbol

x = x(:); % Ensure column vector

% Identify unique symbols and their probabilities
symbols = unique(x);
p = zeros(length(symbols),1);
for i = 1:length(symbols)
    p(i) = sum(x == symbols(i)) / length(x);
end

% Remove zero probabilities
p(p == 0) = [];

% Compute entropy (base-2)
H = -sum(p .* log2(p));

end
