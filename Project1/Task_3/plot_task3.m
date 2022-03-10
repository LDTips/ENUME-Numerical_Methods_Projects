function plot_task3()
    z = task_3(3); % Initialise
    z = task_3_solve(z, 'j'); % Jacobi solve
    jacobi_error = z.iteration_error; % Store our iteration error vector
    [~, n1] = size(z.iteration_error); % Get the amount of iterations for jacobi
    % [columns, rows]. We need only row amount
    % Needed for graph
    clear z; % We reuse z for gauss seidel

    z = task_3(3);
    z = task_3_solve(z, 'g'); % Gauss solve
    gauss_error = z.iteration_error; % Store our iteration error vector
    [~, n2] = size(z.iteration_error); % Get the amount of iterations for gauss
    % [columns, rows]. We need only row amount
    clear z;
    % Plot the graph
    plot(1:n1, jacobi_error, 1:n2, gauss_error);
    legend('jacobi', 'gauss-seidel');
end