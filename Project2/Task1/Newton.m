% It is advised to choose the following intervals for considered function:
% [5.5, 8] ; [8, 10.5]
% Otherwise method will probably diverge
% Highest intervals should be [a, 8] ; [8, b]
function [iter_count, zero_point] = Newton(initial_guess, a, b)
    if (initial_guess < a || initial_guess > b) % The guess must be chosen in the interval
        error("Incorrect initial guess chosen");
    end
    if ( (a < 2 || a > 12) || (b < 2 || b > 12) ) % a and b need to be within [2,12] interval
        error("Chosen interval is not within [2,12] interval");
    end
    precision = 1e-05; % We stop when our zero point is within this tolerance
    zero_point = initial_guess;
    err = inf;
    iter_count = 0;
    % Display initial guess and the chosen interval:
    fprintf("Initial guess: %.5f \t ; chosen interval: [%.3f, %.3f]\n\n", initial_guess, a, b); 
    % Set up for pretty formatting:
    fprintf("Iteration\t\t x value\t\t f(x)\t\t abs(f(zero_point))\n");
    while err > precision
        % As in the formulas from the report
        previous_zero = zero_point; % Store the previous zero point
        zero_point = previous_zero - (eval_f_x(zero_point)/ eval_f_dx(zero_point));
        % xn+1     =    xn         -        f(x)          /         f'(x)
        if (zero_point < a || zero_point > b)
            error('Local divergence detected @ iteration %d ; x = %.5f outside [%.3f, %.3f] interval', iter_count+1, zero_point, a, b);
        end
        err = abs(eval_f_x(zero_point)); % Distance from 0
        iter_count = iter_count + 1;

        % Print the values of the iteration. They are lined up with the
        % above-mentioned fprintf (pretty formatting)
        fprintf("\t%d \t\t %.10f \t %.10f \t %.10f\n", iter_count, zero_point, eval_f_x(zero_point), err);
    end
    format long; % Display the zero point in the 'long' format
end