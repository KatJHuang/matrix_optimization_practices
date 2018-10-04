function [A] = get_freq(J)
    A = bsxfun(@rdivide, J, sum(J, 1));
end