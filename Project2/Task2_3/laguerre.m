function [iter_count, xk] = laguerre(xk) % xk being initial guess, zero point
    format long;
    n = 4; % Order of our polynomial
    iter_count = 0;
    tolerance = 1e-05;
    fprintf("Iteration ; z_min ; xk ; f(xk)\n"); % Define contents of rows
    pretty_print(iter_count, NaN, xk); % Special NaN call to get info for
    % If our initial_guess xk somehow is 100% correct, the
    % 'while' loop does not start. iter_count is returned to be 0
    % We stop iterating when our result is within our tolerance, or
    % divergence occurs (too big iteration count)
    while iter_count < 1000 && abs(eval_f_x(xk)) > tolerance
        % For this algorithm I decided to include denominators separatly
        % For better readability (these are way longer than for MM2)
        denominator_neg = eval_f_dx(xk) - ...
            sqrt( (n-1)*((n-1)*(eval_f_dx(xk)^2) - n*eval_f_x(xk)*eval_f_d2x(xk)) );
        denominator_pos = eval_f_dx(xk) + ...
            sqrt( (n-1)*((n-1)*(eval_f_dx(xk)^2) - n*eval_f_x(xk)*eval_f_d2x(xk)) );

        % Define two different x_k+1 values, as in the formula
        x_neg = (n * eval_f_x(xk)) / denominator_neg;
        x_pos = (n * eval_f_x(xk)) / denominator_pos;

        if abs(denominator_pos) > abs(denominator_neg) % We choose x_k+1 with bigger denominator
            x_min = x_pos;
        else % Otherwise, z_neg has smaller abs value of the denominator
            x_min = x_neg;
        end

        xk = xk - x_min; % Define new zero point for next iteration
        iter_count = iter_count + 1; % bump+ the iteration count
        
        pretty_print(iter_count, x_min, xk); % Prints current iteration info
        % fprintf("%d ; %.10f ; %.10f ; %.10f\n", iter_count, x_min, xk, eval_f_x(xk));
    end

    if (abs(imag(xk)) < 1e-05) % If the imaginary part is negligible, we remove it from our zero point
        xk = real(xk);
    end
    disp('zero_point is '); disp(xk); % Display obtained zero_point
end