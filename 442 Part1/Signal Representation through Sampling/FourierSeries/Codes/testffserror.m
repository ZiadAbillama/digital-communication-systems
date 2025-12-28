% Test: test_ffs_error.m
%Case 1: Fixed T, vary n
T_fixed  = 2;
n_values = [1 10 20 40 80];
E_n      = zeros(size(n_values));

t  = linspace(-T_fixed/2, T_fixed/2, 2000);      % one full period window 
xt = sign(sin(2*pi*2*t));                        % signal on that window

for i = 1:numel(n_values)
    n = n_values(i);
    [xhat, ~] = ffs(xt, t, n, T_fixed);          % integrates over [-T/2, T/2] 
    E_n(i) = trapz(t, abs(xt - xhat).^2);        
end

%Case 2: Fixed n, vary T 
n_fixed = 50;
T_values = [0.6 0.9 1.3 2.5 3.7]; %as T gets further and further away error will increase
E_T = zeros(size(T_values));

for i = 1:numel(T_values)
    T = T_values(i);
    t  = linspace(-T/2, T/2, 2000);              % one full period (varies with T)
    xt = sign(sin(2*pi*2*t));
    [xhat, ~] = ffs(xt, t, n_fixed, T);
    E_T(i) = trapz(t, abs(xt - xhat).^2);
end

figure('Name','Fourier Series Error Analysis','NumberTitle','off');
subplot(1,2,1);
semilogy(n_values, E_n, 'or-','LineWidth',1.5); grid on;
xlabel('n'); ylabel('E(n,T)'); title(sprintf('Error vs n (T=%.1f s)', T_fixed));

subplot(1,2,2);
semilogy(T_values, E_T, 'ob-','LineWidth',1.5); grid on;
xlabel('T (s)'); ylabel('E(n,T)'); title(sprintf('Error vs T (n=%d)', n_fixed));

sgtitle('Finite Fourier Series Approximation Error Analysis');

disp('--- Error vs n (T fixed) ---');
disp(table(n_values.', E_n.', 'VariableNames', {'n','Error_E'}));
disp('--- Error vs T (n fixed) ---');
disp(table(T_values.', E_T.', 'VariableNames', {'T','Error_E'}));
