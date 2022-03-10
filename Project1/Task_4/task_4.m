classdef task_4
    properties
        A; % Matrix we will find eigenvalues from
        tolerance = 1e-06; % Tolerance, per task definition
        n = 5; % Matrix dim, per task definition
        eig_vector; % Eigenvalue matrix
    end
    methods (Access = 'protected')
        % Below methods are the functions found in the book:
        % Piotr Tatjewski - Numerical Methods
        % ----------------------------------------------------------
        % Page 68 - QR factorisation (thin). Full rank matrices only
        % Slightly modified due to the usage of a class for variable storage
        function [obj, Q, R] = QR_factor(obj)
            d = zeros(1, obj.n);
            Q = zeros(obj.n, obj.n);
            R = zeros(obj.n, obj.n);
            % Factorization with orthogonal columns of Q
            for i = 1 : obj.n
                Q(:, i) = obj.A(:, i);
                R(i, i) = 1;
                d(i) = Q(:, i)' * Q(:, i); 
                for j = (i + 1) : obj.n
                    R(i, j) = (Q(:, i)' * obj.A(:, j)) / d(i);
                    obj.A(:, j) = obj.A(:, j) - R(i, j) * Q(:, i);
                end
            end
            % Column normalisation (columns of Q orhonormal)
            for i = 1 : obj.n
                dd = norm(Q(:, i));
                Q(:, i) = Q(:, i) / dd;
                R(i, i:obj.n) = R(i, i:obj.n) * dd;
            end
        end

        % Page 77 - QR eigenvalue algorithm without shifts
        function obj = eig_QR_noShift(obj)
            % tolerance is for '0' elements. Note the display of A matrix
            % After the function completes
            % Tolerance is in the task_4 class, predefined
            i = 1;
            max_iter = 10000;
            while i <= max_iter && max(max(obj.A - diag(diag(obj.A)))) > obj.tolerance
                [obj, Q, R] = QR_factor(obj);
                obj.A = R * Q; % In the end, obj.A will be a matrix with
                % Eigenvalues on the diagonal
                i = i + 1;
            end
            if i > max_iter
                error('Maximum number of iterations @ eig_QR_noShift. Exiting');
            else % To show how many iterations were needed
                disp('Method converged. Iteration count: ');
                disp(i);
            end
            obj.eig_vector = diag(obj.A); % We store the eigenvalues from
            % the diagonal to the eigenvector for better representation
        end

        % Page 77 and 78 - Eigenvalue QR method with shifts
        function obj = eig_QR_shift(obj)
            INITIALsubmat = obj.A;
            max_iter = 10000;
            total_i = 0;
            for k = obj.n:-1 : 2
                DK = INITIALsubmat; % DK - initial matrix
                i = 0;
                while i <= max_iter && max(abs(DK(k, 1:k-1))) > obj.tolerance
                    DD = DK(k-1:k, k-1:k); % Bottom right 2x2 submatrix
                    coeff = [1, -(DD(1, 1) + DD(2, 2)), DD(2, 2) * DD(1, 1) - DD(2 ,1) * DD(1, 2)];
                    eig_roots = roots(coeff); % Solve a quadratic equation for
                    % coeff = [a, b, c]
                    % Choose the eigenvalue closer to bottom right element
                    % of 2x2 submatrix (eigenvalue has smaller distance to)
                    if abs(eig_roots(1) - DD(2, 2)) < abs(eig_roots(2) - DD(2, 2))
                        shift = eig_roots(1);
                    else
                        shift = eig_roots(2);
                    end
                    DP = DK - eye(k) * shift; % Shifted matrix
                    [Q1, R1] = external_QR_factor(DP); % QR decomposition
                    DK = R1 * Q1 + eye(k) * shift;
                    i = i + 1;
                end
                total_i = total_i + i; % Check how many iterations were needed
                if i > max_iter
                    error('Maximum number of iterations @ eig_QR_Shift. Exiting')
                end
                obj.eig_vector(k) = DK(k,k);
                if k > 2
                    INITIALsubmat = DK(1:k-1, 1:k-1);
                else
                    obj.eig_vector(1) = DK(1, 1); % Last eigval
                end
            end
            disp('Total iterations for every submatrix: ');
            disp(total_i);
        end

    end
    methods
        function obj = task_4()
            obj.eig_vector = zeros(obj.n, 1);
            obj.A = randi([-5, 5], obj.n, obj.n); % Generate random matrix 5x5
            % With values from -55 to 55
            obj.A = obj.A + (obj.A)'; % Turn matrix into symmetric

            % I didnt get a single matrix from randi that wasn't full rank
            % But better safe than sorry
            if rank(obj.A) ~= 5
                error('Error: Matrix generated not full rank');
            end
        end
        function obj = eig_QR(obj, doShift)
            % Display matrix A before eigenvalue computation
            disp('Matrix A before finding eigenvals: ');
            disp(obj.A);
            % Do shift:
            if doShift == true
                disp('SHIFT mode');
                obj = eig_QR_shift(obj);
            % Dont shift:
            elseif doShift == false
                disp('NO SHIFT mode');
                obj = eig_QR_noShift(obj);
            end
            disp('Result eigenvector: ');
            disp(obj.eig_vector);
            format long; % Used to show every element apart from diagonals
            % Are smaller than tolerance
            disp('Matrix A after finding eigenvals: ');
            disp(obj.A);
        end
    end
end