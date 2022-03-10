function y = eval_f_dx(x)
    % Interval error handling is done in different functions.
    %if x < 2 || x > 12 % We shouldn't take values outside interval [2,12]
        %error('Tried to evaluate function outside defined interval.');
    %end
    y = (6*cos(x))/5 + 2/(x + 2); % Derivative calculated beforehand using:
    % syms x
    % diff(1.2 * sin(x) + 2 * log(x+2) - 5)
end