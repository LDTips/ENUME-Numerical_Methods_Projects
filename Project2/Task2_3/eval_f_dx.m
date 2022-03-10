function y = eval_f_dx(x)
    format long;
    % a0 = 1;
    a = [2, 5, 5, -2]; % Task description a numeration will be the same
    % Hence a is reversed. 1 is a(1) 2 is a(2), 5 is a(3) etc.
    y = a(4) * 4*x^3 + a(3) * 3*x^2 + a(2) * 2*x + a(1); % Manually computed derivative (1st)
end