function Y = ttm(X, A, n)
%TTM Performs n-mode multiplication of a tensor with a matrix.
%
%   Y = ht.ttm(X, A, n)
%
%   Computes Y = A ×_n X, i.e., multiplies tensor X along its n-th mode
%   by matrix A (from the left). Equivalent to:
%       Y = tensor-times-matrix along mode n
%
%   Inputs:
%       X - Full tensor (numeric multidimensional array)
%       A - Matrix (size: [R, size(X,n)])
%       n - Mode (positive integer, 1-based indexing)
%
%   Output:
%       Y - Resulting tensor after mode-n multiplication

    sz = size(X);
    sz(end+1:n) = 1;
    N = max(ndims(X),n);

    % Permute to bring mode-n to front
    order = [n, 1:n-1, n+1:N];
    X_perm = permute(X, order);
    X_mat = reshape(X_perm, sz(n), []);  % matricize

    % Multiply
    Y_mat = A * X_mat;

    % Reshape back
    new_sz = [size(A,1), sz([1:n-1, n+1:end])];
    Y_perm = reshape(Y_mat, new_sz);

    % Inverse permutation
    inv_order(order) = 1:N;
    Y = permute(Y_perm, inv_order);
end