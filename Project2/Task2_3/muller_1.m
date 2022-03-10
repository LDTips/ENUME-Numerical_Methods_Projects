function [iter_count, x2] = muller_1(x0, x1, x2)
    format long;
    % x2 is our zero point.
    tolerance = 1e-05;
    iter_count = 0;
    fprintf("Iteration ; z_min ; x2 ; f(x2)\n"); % Define contents of rows
    pretty_print(iter_count, NaN, x2); % Special NaN call to get info for
    % 0th iteration before starting the loop

    % If our initial_guess x2 somehow is 100% correct, the
    % 'while' loop does not start. iter_count is returned to be 0
    % We stop iterating when our result is within our tolerance, or
    % divergence occurs (too big iteration count)
    while iter_count < 1000 && abs(eval_f_x(x2)) > tolerance
        % x2 is assumed to be the initial guess
        % Define z0 and z1 per definition:
        z0 = x0 - x2;
        z1 = x1 - x2;

        % Computed a and b values below. c is trivial:
        %a = (-z0 * f(x1) + z0 * f(x2) + z1 * f(x0) - z1 * f(x2)) / (z0^2 * z1 - z0 * z1^2)
        a = (-z0 * eval_f_x(x1) + z0 * eval_f_x(x2) + z1 * eval_f_x(x0) - z1 * eval_f_x(x2)) / ...
            (z0 * z1 * (z0 - z1));
        %b = (z0^2 * (f(x1) - f(x2)) + z1^2 * (f(x2) - f(x0))) / z0*z1 * (z0-z1)
        b = (z0^2 * (eval_f_x(x1) - eval_f_x(x2)) + z1^2 * (eval_f_x(x2) - eval_f_x(x0))) / ...
            (z0 * z1 * (z0 - z1));
        c = eval_f_x(x2);

        % Define z- and x+
        z_neg = -2*c / (b - sqrt(b.^2 - 4*a*c)); % z-
        z_pos = -2*c / (b + sqrt(b.^2 - 4*a*c)); % z+

        if abs(z_neg) > abs(z_pos) % That means z_pos has smaller denominator
            z_min = z_pos;
        else % Otherwise, z_neg has smaller abs value of the denominator
            z_min = z_neg;
        end
        % We need to choose new x0 and x1. The ones being closest to x3
        x3 = x2 + z_min; % Start the iteration with a new zero_point
        if (abs(x0 - x3) > abs(x1 - x3) && abs(x0 - x3) > abs(x2 - x3)) % That means x0 is farthest away
            % So we do not assign it
            x0 = x1;
            x1 = x2;
        elseif (abs(x1 - x3) > abs(x0 - x3) && abs(x1 - x3) > abs(x2 - x3)) % x1 is farthest away
            %x0 = x0;
            x1 = x2;
        end
        iter_count = iter_count + 1; % bump+ the iteration count
        x2 = x3;
        pretty_print(iter_count, z_min, x2); % Prints current iteration info
        %fprintf("%d ; %.10f ; %.10f ; %.10f\n", iter_count, z_min, x2, eval_f_x(x2));
    end

    if (abs(imag(x2)) < 1e-05) % If the imaginary part is negligible, we remove it from our zero point
        x2 = real(x2);
    end
    disp('zero_point is '); disp(x2); % Display obtained zero_point
end   
