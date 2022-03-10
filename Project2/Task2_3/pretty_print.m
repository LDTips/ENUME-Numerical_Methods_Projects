function pretty_print(iter, a_min, x) % Method realises pretty printout of every iteration for MM2, MM1 and laguerre
    if isnan(a_min) % Special call for 0th iteration is indicated by a_min = NaN
        fprintf("0 ; N/A ; %.5f + %.5fi ; %.5f + %.5fi\n", real(x), imag(x), real(eval_f_x(x)), imag(eval_f_x(x)));
    else
        fprintf("%d ; %.5f + %.5fi ; %.5f + %.5fi ; %.5f + %.5fi\n", ...
            iter, real(a_min), imag(a_min), real(x), imag(x), real(eval_f_x(x)), imag(eval_f_x(x)));
    end
end
% This method was written to reduce number of code lines in the algorithm
% That are not directly related to the algorithm itself (the calculations)