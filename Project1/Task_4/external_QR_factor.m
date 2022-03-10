% Page 68 - QR factorisation (thin). Full rank matrices only
 % Slightly modified due to the usage of a class for variable storage
 % Numerical Methods - Piotr Tatjewski

 % This is a modification of already existing function in a class
 % It is need so that QR factorisation can be done for some arbitrary
 % matrix, rather than only the matrix found inside the class
function [Q, R] = external_QR_factor(D)
    [m, n] = size(D);
    d = zeros(1, n);
    Q = zeros(m, n);
    R = zeros(n, n);
    % Factorization with orthogonal columns of Q
    for i = 1 : n
        Q(:, i) = D(:, i);
        R(i, i) = 1;
        d(i) = Q(:, i)' * Q(:, i); 
        for j = (i + 1) : n
            R(i, j) = (Q(:, i)' * D(:, j)) / d(i);
            D(:, j) = D(:, j) - R(i, j) * Q(:, i);
        end
    end
    % Column normalisation (columns of Q orhonormal)
    for i = 1 : n
        dd = norm(Q(:, i));
        Q(:, i) = Q(:, i) / dd;
        R(i, i:n) = R(i, i:n) * dd;
    end
end