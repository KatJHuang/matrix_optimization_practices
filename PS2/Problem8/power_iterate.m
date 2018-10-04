function [u1, error] = power_iterate(A)
    % initialize hypothesis vector x with L2 norm = 1
    x_init = rand(size(A, 1), 1);
    x_init = x_init / norm(x_init);

    max_iter = 10;
    
    x = zeros(size(A, 1), max_iter);
    x(:, 1) = x_init;
    
    y = zeros(size(A, 1), max_iter);
    lambda = zeros(1, max_iter);
    
    error = zeros(1, max_iter);
    
    for ii = 1:max_iter-1
       y(:, ii + 1) = A * x(:, ii);
       x(:, ii + 1) = y(:, ii + 1) / norm(y(:, ii + 1));
       lambda(ii + 1) = x(:, ii + 1)' * A * x(:, ii + 1);
       
       error(:, ii + 1) = norm(A * x(:, ii + 1) -  x(:, ii + 1));
    end
    
    error = log(error(2:max_iter));
    u1 = x(:, max_iter);
end