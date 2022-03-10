function y = eval_f_x(x)
    % Interval error handling is done in different functions.
    %if x < 2 || x > 12 % We shouldn't take values outside interval [2,12]
        %error('Tried to evaluate function outside defined interval');
    %end
    y = 1.2 * sin(x) + 2 * log(x+2) - 5; % As in the task description
end