function plot_im_roots()
    a = roots( [-2 5 5 2 1] );
    plot(real(a), imag(a), 'go');
    grid on;
    xlabel('Re');
    ylabel('Im');
end