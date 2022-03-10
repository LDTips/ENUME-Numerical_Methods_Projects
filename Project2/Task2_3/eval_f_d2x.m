function y = eval_f_d2x(x)
    format long;
    % a0 = 1;
    a = [2, 5, 5, -2]; % Task description a numeration will be the same
    % Hence a is reversed. 1 is a(1) 2 is a(2), 5 is a(3) etc.
    y = a(4) * 12*x^2 + a(3) * 6*x + a(2) * 2; % Manually computed derivative (2nd)
end