classdef task_2
    properties
        A; % Main matrix
        b; % = matrix
        x; % Result matrix. x_1 , x_2, ... , x_n
        res_error; % Result error
        n; % Matrix size. n from the task
    end
    methods (Access = 'private')
        function obj = task_a_create_arr(obj)
            % Fill matrix A:
            [row, col] = size(obj.A); % Used so we can easily address places in the matrix
            for i = 1 : row
                for j = 1 : col
                    % Filling of matrix A in accordance to formula given in
                    % the .pdf
                    if i == j
                        obj.A(i,j) = 9;
                    elseif  (i == j - 1 || i == j + 1)
                        obj.A(i,j) = 3;
                    else 
                        obj.A(i,j) = 0;
                    end
                end
                % Filling matrix b. We dont need row, col addressing
                % Since b is actually a vector
                % Filling in accordance to formula:
                obj.b(i) = 1.5 + 0.5 * i;
            end
        end
        function obj = task_b_create_arr(obj)
            [row, col] = size(obj.A); % Same as task_a_create_arr
            for i = 1 : row
                % Fill matrix A in accordance with task description
                % General formula
                for j = 1 : col
                    obj.A(i,j) = 2/(3*(i + j - 1));
                end
                % Fill vector b in accordance with task description
                % if i is even
                if mod(i, 2) == 0
                    obj.b(i) = 1/i;
                % if i is odd
                else 
                    obj.b(i) = 0;
                end
            end
        end
        
        function obj = gauss_partial_pivot(obj)
            for i = 1 : obj.n
                [local_max, max_row] = max(obj.A(i:obj.n, i)); % Maximum value from column i
                % and rows from i to n
                % we mustn't check rows already processed, that is every
                % row above ith
                max_row = max_row + i - 1;
                % This is done due to how max_row is saved after max() is
                % called
                % e.g. if we check rows 5:n, and the max is in 7th row
                % Then max_row is going to be 3, instead of 7
                % In other words - it gives a row number relative to the
                % chosen starting row rather than relative to 1st row

                % used for pivoting
                % refer to max() function documentation
                if (local_max == 0)
                    disp("Invalid matrix for Gaussian elimination");
                    return
                end
                % Perform pivoting:
                if (max_row ~= i) % If the biggest value
                    % Is not in the current row
                    % Swap rows:
                    obj.A([max_row, i], :) = obj.A([i, max_row], :);
                    obj.b([max_row, i]) = obj.b([i, max_row]);
                    % : as 2nd argument is not needed for b, as it's only a
                    % vector. For A we need to indicate we don't do
                    % anything with columns
                end
                % Perform Gaussian elimination:
                for j = (i + 1) : obj.n
                    l = obj.A(j, i) / obj.A(i, i); % Define curr. row multiplier
                    % In accordance with the definition
                    if l ~= 0 % If l = 0 we can skip these operations
                        % Every value will remain unchanged for l == 0
                        for k = (i + 1) : obj.n
                            obj.A(j, k) = obj.A(j, k) - l * obj.A(i, k);
                            % Change every value in a given row in
                            % accordance with the formula
                        end
                        obj.b(j) = obj.b(j) - l * obj.b(j);
                        % For b we only need to change a single object. 
                        % b is one column
                    end
                end
            end
            % Obtain results from upper-triangular matrix. 
            % Put them into x vector
            obj.x(obj.n) = obj.b(obj.n) / obj.A(obj.n, obj.n); 
            % Trivial x_n formula. From the book
            % Note that we solve the matrix "from the bottom"
            for k = (obj.n - 1):-1 : 1
                summation = 0;
                for j = (k + 1) : obj.n % Sum from j = k+1 up to n
                    summation = summation + obj.A(k,j) * obj.x(j);
                end
                obj.x(k) = (obj.b(k) - summation) / obj.A(k,k);
                % Obtain x_k as written in the formula
            end
            % All above is done according to the formula from the book
            % Subindexes are the same as in the formula (k and j)
        end
    end

    methods
        function obj = task_2(n) % Constructor. Fills matrices with 0, assigns n
            obj.n = n;
            obj.A = zeros(n, n); % n x n matrix
            obj.b = zeros(n, 1); % n x 1 matrix
            obj.x = zeros(n, 1); % n x 1 matrix
            obj.res_error = zeros(n, 1); % n x 1 matrix. 1 value for every result
        end
        function obj = task_2a(obj)
            obj = task_a_create_arr(obj); % Create the array for the task a
            obj = gauss_partial_pivot(obj); % Do the gaussian elimination
            obj = calc_err(obj); % Calculate the result error
        end
        function obj = task_2b(obj)
            obj = task_b_create_arr(obj); % Create the array for the task b
            obj = gauss_partial_pivot(obj); % Do the gaussian elimination
            obj = calc_err(obj); % Calculate the result error
        end
        function display_results(obj) % Purely for display
            % We dont need to store these in any new variables
            format compact
            obj.x
            obj.res_error
        end
        function obj = iterative_refinement(obj)
            % All done according to formulas from book
            % First, we need to solve a different system of equations
            % Instead of Ax = b, we solve A(delta x) = r where r is
            % residuum, delta x is error. Residuum is res_error
            correction_obj = obj; % Copy current object
            correction_obj.b = correction_obj.res_error; % Switch b with delta x
            correction_obj = gauss_partial_pivot(correction_obj); % Solve
            for i = 1 : obj.n
                obj.x(i) = obj.x(i) - correction_obj.x(i);
                % x_2 = x_1 - delta x
                % We override our obj with the new
                % iteration. We dont need old iteration
            end
            obj = calc_err(obj); % Recalc the error
        end
        function obj = calc_err(obj)
            % We use the formula: r = Ax - b
            for i = 1 : obj.n
                summation = 0; % Assign temp variable
                for j = 1 : obj.n
                    summation = summation + obj.A(i,j) * obj.x(j); % Ax
                    % We need iterations so that we include every value in
                    % a given row. i (row) is const. in this local for loop
                end
                obj.res_error(i) = summation - obj.b(i); % Ax - b
            end
        end
    end
end