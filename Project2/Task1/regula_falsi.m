% It is advised to choose the following intervals for considered function:
% [2, 8] ; [8, 12]
% It will also work for [2, 12], but it will only find one zero point
function [iter_count, c] = regula_falsi(a, b) % a,b - interval to operate on. c is our zero point
    if ( (a < 2 || a > 12) || (b < 2 || b > 12)) % a and b need to be within [2,12] interval
        error("Chosen interval is not within [2,12] interval");
    end
    if ( eval_f_x(a) * eval_f_x(b) >= 0) % Otherwise the method might diverge
        % We are unsure if there are any zero points in the interval if
        % above is true
        error("Function has no different signs at interval boundaries. f(a) * f(b) >= 0");
    end
    fprintf("Chosen interval: [%.3f, %.3f]\n\n", a, b);
    fprintf("Iteration\t\t x value\t\t f(x)\t\t\t abs(f(x))\t\t interval width abs(a-b)\n");
    iter_count = 0;
    err = inf;
    accuracy = 1e-05;
    while err > accuracy
        % we will use the formula from the book:
        % Calculate c point:
        c = (a * eval_f_x(b) - b * eval_f_x(a)) / (eval_f_x(b) - eval_f_x(a));
       %c = (a *     f(b)    - b *    f(a)    ) / (   f(b)     -    f(a)    )

        % Define new interval for next iteration. Per book formula:
        if (eval_f_x(a) * eval_f_x(c) < 0)
            % a = a. unchanged
            b = c;
        elseif (eval_f_x(c) * eval_f_x(b) < 0)
            a = c;
            % b = b. unchanged
        end
        
        err = abs(eval_f_x(c)); % Distance from 0
        iter_count = iter_count + 1;

        fprintf("\t%d \t\t %.10f \t %.10f \t\t %.10f \t\t %.10f\n", iter_count, c, eval_f_x(c), err, abs(a-b));
    end
    format long; % Display the zero point in proper format
end