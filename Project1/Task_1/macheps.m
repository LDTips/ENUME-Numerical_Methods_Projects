macheps_ = 1; % Define the variable. We divide it by /2 at the first condition checking, hence 1.00

while (1.00 + (macheps_ / 2) > 1.00) % From the definition of epsilon. 1 + eps = 1
    macheps_ = macheps_ / 2; % We divide macheps by 2 to search for smaller number
    % This operation is equivalent to shifting mantissa to the right
    % 2^-1 -> 2^-2 -> ...
end
% macheps is named macheps_, because variable names can't be the same as
% script names...