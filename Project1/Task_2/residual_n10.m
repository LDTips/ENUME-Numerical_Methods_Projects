function residual_n10(task_letter, correction_count)
    z = task_2(10); % Predefined, we know n = 10
    if task_letter == 'a'
        z = task_2a(z);
    elseif task_letter == 'b'
        z = task_2b(z);
    end

    % Display results before residual correction(s)
    disp('Results before all corrections:');
    for i = 1 : z.n
        fprintf("x_%d= \t %.25f\n", i, z.x(i)); %x(i) to 25th precision
        % Prints current iteration number and corresponding x vector value
    end

    % residual corrections:
    % Prints max(abs(error)) before any iterative refinement
    disp('Correction information:');
    fprintf("correction \t max(abs(error))\n");
    fprintf("\t0 \t %.25f\n", max(abs(z.res_error)));
    for i = 1 : correction_count
        z = iterative_refinement(z); % Correction
        fprintf("\t%d \t %.25f\n", i, max(abs(z.res_error)));
        % Prints current iteration number and maximum absolute value from
        % res_error vector. z is changed every iteration
    end

    % Prints the results after all corrections were done
    disp('after all corrections:');
    for i = 1 : z.n
        fprintf("x_%d= \t %.25f\n", i, z.x(i));
    end

    % Prints 1-condition and 2-condition number for matrix A.
    % Used to show that matrix A from subtask b is ill-conditioned
    disp('1-condition:')
    disp(cond(z.A, 1));
    disp('2-condition:');
    disp(cond(z.A, 2));
    
end