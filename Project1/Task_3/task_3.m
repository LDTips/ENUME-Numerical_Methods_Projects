classdef task_3
    % This class has elements copied from task_2 class. Namely:
    % Variables; task_b and task_a array creation; error calculation
    properties
        A; % Main matrix
        b; % = matrix
        x; % Result matrix. x_1 , x_2, ... , x_n
        res_error; % Result error
        n; % matrix size

        accuracy; % Accuracy defined as 10^-10
        iteration_error; % Table of the errors every iteration
    end
    
    methods (Access = 'private')
        function obj = jacobi(obj)
            obj.iteration_error(1) = inf; % So that the while loop starts
            % The values will become legitimate after first while loop
            limit = 10000;
            k = 1;
            % error = inf;
            % obj is x(k)
            while k <= limit && obj.iteration_error(k) > obj.accuracy
                % We end either when we:
                % reach the iteration limit, or the error is satisfied
                x_prev = obj.x;         
                for i = 1 : obj.n
                    summation = 0;
                    for j = 1 : obj.n
                        if j ~= i % Summation requirement
                            summation = summation + (obj.A(i,j) * x_prev(j));
                            % Do the summation as in the formula
                        end
                    end
                    obj.x(i) = (obj.b(i) - summation) / obj.A(i,i);
                    % Plug in the rest of variables as it is in the formula
                end
                k = k + 1;
                obj = calc_err(obj); % Recalc r (error). Used next line
                obj.iteration_error(k) = norm(obj.res_error); % STOP TEST formula
            end
            obj.iteration_error = obj.iteration_error(2:end);
            % remove our dummy inf at 1st position
            % Display information about the performed algorithm
            disp('Iteration count: ');
            disp(k-1);
            format long;
            disp('Stopped at error: ');
            disp(obj.iteration_error(end));
            disp('Result matrix: ')
            disp(obj.x);
        end
        function obj = gauss_seidel(obj)
           obj.iteration_error(1) = inf; % So that the while loop starts
            % The values will become legitimate after first while loop
            limit = 10000;
            k = 1;
            % obj is x(k)
            while k <= limit && obj.iteration_error(k) > obj.accuracy
                % We end either when we:
                % reach the iteration limit, or the error is satisfied
                x_prev = obj.x;         
                for i = 1 : obj.n
                    summation = 0;
                    for j = 1 : i - 1 % Bounds in the formula
                        summation = summation + (obj.A(i,j) * obj.x(j));
                        % As in the formula
                    end
                    for j = i + 1 : obj.n % Bounds in the formula
                        summation = summation + (obj.A(i,j) * x_prev(j));
                        % As in the formula
                    end
                    obj.x(i) = (obj.b(i) - summation) / obj.A(i,i);
                    % Plug in the rest of variables as it is in the formula
                end
                k = k + 1;
                obj = calc_err(obj); % Recalc r (error). Used next line
                obj.iteration_error(k) = norm(obj.res_error); % STOP TEST
            end
            obj.iteration_error = obj.iteration_error(2:end);
            % remove our dummy inf at 1st position
            % Display information about the performed algorithm
            disp('Iteration count: ');
            disp(k-1);
            format long;
            disp('Stopped at error: ');
            disp(obj.iteration_error(end));
            disp('Result matrix: ')
            disp(obj.x);
        end
        function obj = task_2a_create_arr(obj) % Method from task_2 class
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
        function obj = task_2b_create_arr(obj) % Method from task_2 class
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
        function obj = task_3_create_arr(obj)
            obj.A = [6 2 1 -1; 4 -10 2 -1; 2 -1 8 -1; 5 -2 1 -8]; % predefined
            obj.b = [6; 8; 0; 2]; % Predefined
        end
    end

    methods
        function obj = task_3(task_number) % Method from task_2 class, changed
            % Constructor. Fills matrices with 0, assigns n
            if task_number == 3
                obj.n = 4; % For task 3, we have 4 by 4 A matrix
            elseif task_number == 2
                obj.n = 10; % For task 2 we have 10 by 10 A matrix
            else 
                return
            end
            obj.A = zeros(obj.n, obj.n); % n x n matrix
            obj.b = zeros(obj.n, 1); % n x 1 matrix
            obj.x = zeros(obj.n, 1); % n x 1 matrix
            obj.res_error = zeros(obj.n, 1); % n x 1 matrix. 1 value for every result

            obj.accuracy = 1e-10;
            obj.iteration_error = [];
        end
        function obj = task_3_solve(obj, type)
            obj = task_3_create_arr(obj);
            if (type == 'j')
                obj = jacobi(obj);
            elseif (type == 'g')
                obj = gauss_seidel(obj);
            end
        end
        function obj = task_2a_solve(obj, type) % Method from task_2 class, changed
            obj = task_2a_create_arr(obj); % Create the array for the task a
            if (type == 'j')
                obj = jacobi(obj);
            elseif (type == 'g')
                obj = gauss_seidel(obj);
            end
        end
        function obj = task_2b_solve(obj, type) % Method from task_2 class, changed
            obj = task_2b_create_arr(obj); % Create the array for the task b
            if (type == 'j')
                obj = jacobi(obj);
            elseif (type == 'g')
                obj = gauss_seidel(obj);
            end
        end
        function obj = calc_err(obj) % Method from task_2 class, changed
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
