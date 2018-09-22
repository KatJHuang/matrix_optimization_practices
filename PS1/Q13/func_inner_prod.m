function [inner_prod] = func_inner_prod(f, g, dx)
%FUNC_INNER_PROD approximates the inner product between two functions
%   Detailed explanation goes here

    inner_prod = (f' * g) * dx;
end

