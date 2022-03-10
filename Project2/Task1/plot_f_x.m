function plot_f_x(step)
    boundaries = [2, 12];
    interval = abs(boundaries(1) - boundaries(2));

    % 'step' is inverted due to the fact that for e.g. step = 0.01
    % We  will need to iterate 100 times for every single whole digit
    % In other words, we will need to iterate 1/step times for every [2,3]
    % [3,4] etc.
    f_x_arr = zeros(1, interval * 1/step); % Preallocate the array.
    for i = 0 : size(f_x_arr, 2) % Iterate from 0 until the f_x_arr is filled
        f_x_arr(i+1) = eval_f_x(2 + i * step); % start from f(2)
        % Then next value will be f(2 + current_iteration * step)
    end
    x_range = 2:step:12; % Define the x_axis values. array of 2 to 12 with
    % difference between every element being 'step'
    mark = (abs(f_x_arr) < 1e-03); %  Define the condition for 0 point marking
    % We cannot use == 0, because our data is generated the way that 0
    % might not be present in the plot at any x considered by this algo
    plot(x_range, f_x_arr, 'k'); % Plot a dark line of f(x)
    % x axis is x_range, f_x_arr are corresponding values
    hold on; % Used for adding the zero points
    plot(x_range(mark), f_x_arr(mark), '*g'); % Add points that satisfy
    % criteria defined by 'mark'. *g == mark points with a green star

    % Added from console after the plot is created:
    % hold on
    % yline(0)
end