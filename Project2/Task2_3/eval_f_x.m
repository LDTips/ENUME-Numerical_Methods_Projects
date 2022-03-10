function y = eval_f_x(x)
    format long;
    a0 = 1;
    a = [2, 5, 5, -2]; % Task description a numeration will be the same
    % Hence a is reversed. 1 is a(1) 2 is a(2), 5 is a(3) etc.
    y = a(4) * x^4 + a(3) * x^3 + a(2) * x^2 + a(1) * x + a0; % As in the task description
    % a0 is not in array, because in matlab arrays start at 1
end