classdef task_2
    properties
        odefun; % Our set of ode equations per task def as anonymous function
        t_interval; % Interval per task def
        x_init; % Initial values per task def
    end

    methods
        function obj = task_2()
            % Values are taken from the task .pdf
            % Anonymous function type:
            obj.odefun = @(t, x) [x(2) + x(1) * (0.5 - x(1)^2 - x(2)^2); -x(1) + x(2) * (0.5 - x(1)^2 - x(2)^2)]; % Explained in report
            obj.t_interval = [0, 15];
            obj.x_init = [0, -0.3]; % x1(0) = 0 ; x2(0) = -0.3 from task desc
        end

        function [t, x] = RK4(obj, h, special_mode)
            t = obj.t_interval(1):h:obj.t_interval(2); % Define our t axis values
            if (exist('special_mode', 'var') == 0) % if special_mode was not passed e.g. in manual fun call
                % Set it to false
                special_mode = false;
            end
            if (special_mode == true) % 5 first values for Adams only
                t = t(1):h:t(5); % Reduce out t array to account for above
            end
            x(:, 1) = obj.x_init; % Assign initial x values
            for n = 1 : (length(t) - 1)
                % All below are formulas from the report for k1, k2, k3, k4 and
                % next iteration y_n+1
                k1 = obj.odefun(t(n), x(:, n));
                k2 = obj.odefun(t(n) + h/2, x(:, n) + h/2 * k1);
                k3 = obj.odefun(t(n) + h/2, x(:, n) + h/2 * k2);
                k4 = obj.odefun(t(n+1), x(:, n) + h * k3); % note: t(n+1) = t + h

                x(:, n+1) = x(:, n) + h * (k1 + 2*k2 + 2*k3 + k4)/6;
            end
        end

        function [t, x] = Adams_PECE_5(obj, h)
            t = obj.t_interval(1):h:obj.t_interval(2); % t axis values
            B_exp = [1901, -2774, 2616, -1274, 251] / 720; % Beta_explicit consts.
            B_imp = [475, 1427, -798, 482, -173, 27] / 1440; % Beta_implicit consts.
            % Note that beta_implicit starts at b*_0, but b_explicit at b_1
            % "Special starting procedure"
            % Explained in the theoretical part of the report
            % It will fetch first five initial values we will use for first step of
            % Adams PECE method
            [~, x] = RK4(obj, h, true); % last arg indicates a special mode for adams pc. Refer to rk4 comments
            for i = 5 : (length(t) - 1)
                sum_temp = 0;
                % Predictor formula. Provided in the report
                for j = 0 : 4
                    % Indexes used can be explained by looking at the "no for" loop
                    % formula provided few lines below
                    sum_temp = sum_temp + B_exp(j + 1) * obj.odefun(t(i - j), x(:, i - j));
                end
                prediction = x(:, i) + h * sum_temp; % Plug in the sum and finalise
                % Without the use of for loop, the full formula would be:
                % prediction = x(:, i) + ...
                                %h * (B_exp(1) * obj.odefun(t(i), x(:, i))...
                                %+ B_exp(2) * obj.odefun(t(i-1), x(:, i-1))...
                                %+ B_exp(3) * obj.odefun(t(i-2), x(:, i-2))...
                                %+ B_exp(4) * obj.odefun(t(i-3), x(:, i-3))...
                                %+ B_exp(5) * obj.odefun(t(i-4), x(:, i-4)));

                % Corrector formula. Provided in the report. First calculation is
                % done outside for loop as it uses prediction variable from before.
                sum_temp = B_imp(1) * obj.odefun(t(i+1), prediction);
                for j = 0 : 4
                    % Indexes used can be explained by looking at the "no for" loop
                    % formula provided few lines below
                    sum_temp = sum_temp + B_imp(j + 2) * obj.odefun(t(i - j), x(:, i - j));
                end
                x(:, i+1) = x(:, i) + h * sum_temp; % Plug in the sum and finalise

                % Without the use of for loop, the full formula would be:
                %t(:, i+1) = t(:, i) + ...
                            %h * B_imp(1) * obj.odefun(t(i+1), prediction) + ...
                            %h * (B_imp(2) * obj.odefun(t(i), x(:, i))...
                                %+ B_imp(3) * obj.odefun(t(i-1), x(:, i-1))...
                                %+ B_imp(4) * obj.odefun(t(i-2), x(:, i-2))...
                                %+ B_imp(5) * obj.odefun(t(i-3), x(:, i-3))...
                                %+ B_imp(6) * obj.odefun(t(i-4), x(:, i-4)));
            end
        end
        function plot_single_alg(obj, h, alg_type, extra_plot)
            % Check the alg_type argument that indicates which alg to use
            if (alg_type == 1)
                [t, x] = RK4(obj, h, false); % alg_type 1 means RK4
                alg_name = "RK4";
            elseif (alg_type == 2)
                [t, x] = Adams_PECE_5(obj, h); % alg_type 2 means PC method
                alg_name = "Adams PC";
            else
                error("Invalid alg_type! Choose 1 for RK4, 2 for Adams PECE.");
            end

            plot(t,x); % Plot x(t)
            legend("x1(t)", "x2(t)"); % Label the lines
            title(alg_name + " x1(t) x2(t) for step = " + h); % Displays algorithm type, graph type and step chosen
            % Same as above, but we plot x1(x2)
            if (extra_plot == true) % This is done so that when we compare functions for different steps
                % In some cases, we might not want to show x1(x2) plot
                figure;
                plot(x(1, :), x(2, :));
                legend("x1(x2)");
                title(alg_name + " x1(x2) for step = " + h); % Displays algorithm type, graph type and step chosen
            end
        end
        function compare_plots(obj, h, iter_count)
            RK4_plot = [figure; figure]; % RK4_plot(1) will be for x1, (2) for 2
            Adams_plot = [figure; figure]; % Same logic here
            % Set plot titles
            figure(RK4_plot(1)); title("RK4 x1(t)");
            figure(RK4_plot(2)); title("RK4 x2(t)");
            figure(Adams_plot(1)); title("Adams x1(t)");
            figure(Adams_plot(2)); title("Adams x2(t)");
            legend_text = strings(1,iter_count); % Used to label every line
            % Show which line corresponds to which step
            for i = 1 : iter_count
                % RK4 plots part:
                [t, x] = RK4(obj, h, false); % Fetch results from RK4 alg
                figure(RK4_plot(1)); % Choose current figure as the 1st RK4
                hold on;
                plot(t, x(1, :)); % Add x1(t) for current h to the plot
                hold off;
                figure(RK4_plot(2));
                hold on;
                plot(t, x(2, :)); % Then add x2(t) for current h to the 2nd plot
                hold off;
                
                % Adams PC plots part:
                [t, x] = RK4(obj, h, false); % Fetch results from RK4 alg
                figure(Adams_plot(1)); % Choose current figure as the 1st RK4
                hold on;
                plot(t, x(1, :)); % Add x1(t) for current h to the plot
                hold off;
                figure(Adams_plot(2));
                hold on;
                plot(t, x(2, :)); % Then add x2(t) for current h to the 2nd plot
                hold off;

                legend_text(i) = "Step = " + h; % Array of strings that will be used to display the legend
                h = h/2; % Decrease the precision by half
            end
            % Add the legends to the plots to show which line corresponds
            % to which step
            figure(RK4_plot(1)); legend(legend_text);
            figure(RK4_plot(2)); legend(legend_text);           
            figure(Adams_plot(1)); legend(legend_text);
            figure(Adams_plot(2)); legend(legend_text);
        end
        % After manually finding the best step, this function is called
        function x1x2_plots(obj, h1, h2)
            %x1(x2) plot part:
            x1x2_plot = [figure; figure]; % 1st figure is for RK4, 2nd for Adams
            % RK4
            figure(x1x2_plot(1)); % Choose the 1st figure (for RK4)
            hold on;
            % Write the title to show initial and final h, as well as
            % algorithm type:
            title("x1(x2) RK4 h1 = " + h1 + " h2 = " + h2);
            [~, x] = RK4(obj, h1, false);
            plot(x(1, :), x(2, :)); % Plot RK4 for initial_h
            [~, x] = RK4(obj, h2, false);
            plot(x(1, :), x(2, :));% Plot RK4 for final_h
            legend("initial-h", "final-h"); % Add the legend
            hold off;
            % Adams PC
            figure(x1x2_plot(2)); % Choose the 2nd figure (for Adams)
            hold on;
            % Write the title to show initial and final h, as well as
            % algorithm type:
            title("x1(x2) Adams PC h1 = " + h1 + " h2 = " + h2); 
            [~, x] = Adams_PECE_5(obj, h1);
            plot(x(1, :), x(2, :)); % Plot Adams for initial_h
            [~, x] = Adams_PECE_5(obj, h2);
            plot(x(1, :), x(2, :)); % Plot Adams for final_h
            legend("h1", "h2"); % Add the legend
            hold off;
        end
        
        % Function that compares initial_h plots vs optimal_h plots
        % Initial_h is h that was one iteration lower than found optimal_h
        function initial_vs_optimal(obj, initial_h, optimal_h)
            RK4_plots = figure;
            Adams_plots = figure;
            x1x2_plots = [figure; figure];
            % Below will be used in title() and legend():
            h_info_string = ["initial-h = " + initial_h, "optimal-h = " + optimal_h];
            % Fetch results of the algorithms:
            % Since t has different sizes for initial and optimal_h, we can
            % not use an array. Same goes for x
            [t1, x1(:, :, 1)] = RK4(obj, initial_h, false);
            [t2, x2(:, :, 1)] = RK4(obj, optimal_h, false);
            [~, x1(:, :, 2)] = Adams_PECE_5(obj, initial_h);
            [~, x2(:, :, 2)] = Adams_PECE_5(obj, optimal_h);
            % We dont need to store t(3) and t(4) because the values will
            % be the same as t axis values from RK4 algorithm
            % Last index of x1 and x2 indicates whether data is for RK4 or
            % Adams. 1 means RK4, 2 means Adams
            % x1 indicates result size for initial_h, x2 indicates result
            % size for optimal_h

            figure(RK4_plots); % RK4 plot operations
            hold on;
            plot(t1, x1(:, :, 1)); % Plot RK4 for initial_h
            plot(t2, x2(:, :, 1)); % Plot RK4 for optimal_h
            % Set title and legend:
            title("x1(t) x2(t) RK4 comparison. " + h_info_string(1) + " " + h_info_string(2));
            legend("x1(t) " + h_info_string(1), "x2(t) " + h_info_string(1), "x1(t) " + h_info_string(2), "x2(t) " + h_info_string(2));
            hold off;

            figure(Adams_plots); % Adams PECE plot operations, analogical to RK4 plot operations
            hold on;
            % Remember that t values are the same as for RK4, so we reuse
            plot(t1, x1(:, :, 2)); % Plot Adams PC for initial_h
            plot(t2, x2(:, :, 2)); % Plot Adams PC for optimal_h
            % Set title and legend:
            title("x1(t) x2(t) Adams PC comparison. " + h_info_string(1) + " " + h_info_string(2));
            legend("x1(t) " + h_info_string(1), "x2(t) " + h_info_string(1), "x1(t) " + h_info_string(2), "x2(t) " + h_info_string(2));
            hold off;
            
            figure(x1x2_plots(1)); % x1x2 plots comparison - RK4 part
            hold on;
            % RK4 needs to have last index = 1. Rest is the same as typical
            % x1(x2) plot. x1 with last index = 1 is RK4, initial_h
            % x2 with last index = 1 is RK4 for optimal_h
            plot(x1(1, :, 1), x1(2, :, 1));
            plot(x2(1, :, 1), x2(2, :, 1));
            % Set title and legend:
            title("x1(x2) RK4 comparison. " + h_info_string(1) + " " + h_info_string(2));
            legend(h_info_string(1), h_info_string(2));
            hold off;

            figure(x1x2_plots(2)); % x1x2 plots comparison - Adams PC part
            hold on;
            % Adams PC needs to have last index = 2. Rest is the same as typical
            % x1(x2) plot. x1 with last index = 2 is Adams PC, initial_h
            % x2 with last index = 2 is Adams PC for optimal_h
            plot(x1(1, :, 2), x1(2, :, 2));
            plot(x2(1, :, 2), x2(2, :, 2));
            % Set title and legend:
            title("x1(x2) Adams PC comparison. " + h_info_string(1) + " " + h_info_string(2));
            legend(h_info_string(1), h_info_string(2));
            hold off;
        end
    end
end