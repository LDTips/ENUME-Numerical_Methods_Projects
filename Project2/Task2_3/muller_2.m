function [iter_count, zero_point] = muller_2(zero_point)
    % zero_point is denoted as x0 in later iterations
    format long; % Suprisingly important instruction. Otherwise the precision will be bad
    tolerance = 1e-05;
    iter_count = 0;
    fprintf("Iteration ; z_min ; x0 ; f(x0)\n"); % Define contents of rows
    pretty_print(iter_count, NaN, zero_point); % Special NaN call for 0th iteration
    % If our initial_guess zero_point somehow is 100% correct, the
    % 'while' loop does not start. iter_count is returned to be 0
    % We stop iterating when our result is within our tolerance, or
    % divergence occurs (too big iteration count)
    while iter_count < 1000 && abs(eval_f_x(zero_point)) > tolerance
        % As in the formula from the book, define a,b,c:
        c = eval_f_x(zero_point);
        b = eval_f_dx(zero_point);
        a = eval_f_d2x(zero_point);

        % Define z- and x+
        z_neg = -2*c / (b - sqrt(b^2 - 2*a*c)); % z-
        z_pos = -2*c / (b + sqrt(b^2 - 2*a*c)); % z+

        % We want the z value with bigger absolute DENOMINATOR value
        if abs(z_neg) > abs(z_pos) % That means z_pos has smaller denominator
            z_min = z_pos;
        else % Otherwise, z_neg has smaller abs value of the denominator
            z_min = z_neg;
        end

        zero_point = zero_point + z_min; % Otherwise start the iteration with a new zero_point
        iter_count = iter_count + 1; % bump+ the iteration count
        
        % Iteration information printout
        pretty_print(iter_count, z_min, zero_point); % Prints current iteration info
        %helper = eval_f_x(zero_point);
        %fprintf("%d ; %.10f ; %.10f ; %.10f + %.10fi\n", iter_count, z_min, zero_point, real(helper), imag(helper));
    end

    if (abs(imag(zero_point)) < 1e-05) % If the imaginary part is negligible, we remove it from our zero point
        zero_point = real(zero_point);
    end
    disp('zero_point is '); disp(zero_point); % Display obtained zero_point
end