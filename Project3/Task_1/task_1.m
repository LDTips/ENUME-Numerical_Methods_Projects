classdef task_1
    properties 
        x; % Experimental x values
        y; % Experimental y values
        a; % Result matrix a (a coefficients)
    end
    % function [Q, R] = external_QR_factor(D) will be used, but is not
    % inside the class

    methods
        % Class initialiser
        function obj = task_1()
            % Assign neccessary x and y values, per task description
            obj.x = -5:5; % x = [-5, -4, -3, [...] , 4, 5];
            obj.y = [-4.9606; -3.3804; -1.4699; -1.1666; 0.4236; 0.1029; ...
                    -0.5303; -4.0483; -11.0280; -21.1417; -33.9458];
            % Transpose x (due to how -5:5 works)
            obj.x = obj.x';
        end

        % This function runs the Least Square Approximation Algorithm with
        % the use of QR factorisation. Result is stored in obj.a
        function [obj, cond_gram, error_norm] = Least_Square_Approx_QR(obj, deg)
            % Note: Before the last line of this algorithm, obj.a is not
            % our 'a' matrix, but 'A' matrix instead.
            % I chose to reuse obj.a as 'A' before storing the actual
            % result to optimise memory use
            obj.a = zeros(size(obj.x, 1), deg+1); % Define empty A
            % Of the size being (length of x) x (degree + 1)
            % Degree 0 is valid, hence +1
            
            % Refer to the computed a matrix at page 5 of the report:
            for i = 1:size(obj.x, 1) % Amount of rows is the amount of data points
                for j = 1:deg+1 % Amount of columns is the degree+1 (because e.g. degree 0 has 1 col)
                    obj.a(i, j) = obj.x(i)^(j-1); % Since arrays in MATLAB start at 1
                    % But we need x^0 power to be considered as well, we use j-1 as exponent
                end
            end
            % Condition number of our Gram's matrix needs to be calculated now as Gram's matrix is discarded
            cond_gram = cond(obj.a' * obj.a); % 

            [Q, R] = external_QR_factor(obj.a); % Perform QR factorisation with the use of thin algorithm
            obj.a = R\(Q'*obj.y); % Ra = Q'y (page 5 in report), solve, store result in a.

            % Error calculation:
            a_poly = poly_a(obj); % Fetch modified a so that polyval can process our 'a' matrix
            computed_values = polyval(a_poly, obj.x); % Evaluate our obtained y values.
            % It evaluates polynomial values with coefficients A_poly at
            % the points from obj.x
            error_norm = norm(computed_values - obj.y); % Residuum, norm(Ax - y)
        end
        
        % This method plots a single line using stored A matrix
        % It can be chosen whether to mark our 'experimental samples'
        function plot_single(obj, markers)
            if (isempty(obj.a) == true) % Least_Square_Approx_QR should be run prior to running of this function
                error("'a' matrix is empty! Unable to plot nonexistent data!");
            end
            a_poly = poly_a(obj); % Obtain a modified A matrix that will be parsed in poly MATLAB functions
            % Our computed a matrix is not fit for the syntax of polyval
            x_axis = -5:.1:5; % Define the x axis values. -5 to 5 with intervals of 0.1
            figure; % Creates an empty figure. 
            %This is done so that there can be many windows with different plots
            plot(x_axis, polyval(a_poly, x_axis)); % Plots our defined x_axis vs values of the computed polynomial
            title("Degree " + (size(obj.a, 1) - 1) ); % Puts a title on the plot to show what degree the plot shows
            grid on;
            if (markers == true) % Passed as argument, whether points should be marked
                hold on;
                plot(obj.x, obj.y, 'LineStyle', 'none', 'Marker', 'o', 'Color', 'red');
                hold off;
            end
        end
        
        % Function below calls all neccessary methods that the task
        % description mentions: Least Square Approximation, error_norm
        % calculation and then plots the computed graph with marked points 
        function obj = all_in_one(obj, deg, display_a)
            [obj, cond_gram, error_norm] = Least_Square_Approx_QR(obj, deg); % Execute least Square approx algorithm
            % Also save the condition number of the Gram's matrix
            % Display all the information:
            fprintf("Degree ; Condition ; Error_norm (csv format)\n%d ; %.5f ; %.10f\n", deg, cond_gram, error_norm); 
            if (display_a == true) % We choose whether to display 'a' matrix or not
                disp("a transposed: [a0, a1, ...]^T (csv format)"); 
                for i = 1:size(obj.a, 1)                 
                    fprintf("%.5f ; ", obj.a(i)'); % Display all matrix a values in a csv format
                end
                fprintf("\n\n");
            end            
            plot_single(obj, true); % Display the plot
        end

        % This function allows to execute multiple least square algorithms
        % for different degress, as well as print neccessary data after
        % every algorithm execution for a degree
        function obj = execute_multiple(obj, min_deg, max_deg, display_a)
            for i = min_deg:max_deg % We execute in the range of the provided min and max deg
                obj = all_in_one(obj, i, display_a); % Refer to the function comments
                pause(0.5); % To avoid fast window pop ups
            end
        end
        % This method returns a new A matrix, one that can be used in poly MATLAB functions
        function a_poly = poly_a(obj) % polyval takes the 'a' coefficients in a reverse way
            a_poly = obj.a;
            a_poly = flip(a_poly); % Reverse the order
            a_poly = a_poly'; % Transpose
        end
    end
end