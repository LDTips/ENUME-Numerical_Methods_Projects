function plot_errors(n, task_letter)
    n_vector = zeros(n, 1); % Vector of N=10,20,... used for plotting later on
    error_vector0 = zeros(n, 1); % Vector with errors for
    %error_vector contains norms of res_error vector for different
    %N, as written in the task description
    for i = 1 : n
        n_vector(i) = 10 * 2^(i-1); % N = 10,20,40,80 = (10 * 2^0), (10* 2^1) ...
        z = task_2(n_vector(i)); % Create empty task_2 object of
        % Size n. n differs every iteration, as it can be seen in
        % line 7
        if task_letter == 'a'
            z = task_2a(z);
        elseif task_letter == 'b'
            z = task_2b(z);
        end
        error_vector0(i) = norm(z.res_error); % Norm of error vector
        % i.e. norm of r = Ax - b
    end
    plot(n_vector, error_vector0);
end