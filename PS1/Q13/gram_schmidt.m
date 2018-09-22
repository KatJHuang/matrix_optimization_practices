function [U] = gram_schmidt(V, dx)
% The vectors v1, ..., vk (columns of matrix V, so that V(:,j) is the jth 
% vector) are replaced by orthonormal vectors (columns of U) which span the
% same subspace
    n = size(V,1); % dimension of each vector
    k = size(V,2); % number of vectors
    
    U = zeros(n,k);
    
    U(:,1) = V(:,1)/sqrt(func_inner_prod(V(:,1), V(:,1), dx)); % normalize the first vector
    for i = 2:k
      U(:,i) = V(:,i);
      for j = 1:i-1
        U(:,i) = U(:,i) - func_inner_prod(U(:,i), U(:,j), dx)...
            / func_inner_prod(U(:,j), U(:,j), dx)...
            * U(:,j);
        
      end
      % normalize after projection onto other bases
      U(:,i) = U(:,i)/sqrt(func_inner_prod(U(:,i), U(:,i), dx)); 
    end
end

