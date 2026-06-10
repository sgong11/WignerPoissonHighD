function nrm = norm(htd, mode)
%NORM Computes the Frobenius norm of an HT tensor.
%
%   nrm = ht.norm(htd)
%   nrm = ht.norm(htd, 'orth')
%
%   Input:
%       htd  - HTD structure with fields .U, .B, .tree, .rank
%       mode - (optional) 'orth' if htd is orthogonalized
%
%   Output:
%       nrm  - Frobenius norm ‖htd‖_F

    if nargin < 2 || ~strcmp(mode, 'orth')
        % General case using inner product
        nrm = sqrt(ht.inner_product(htd, htd));
    else
        % Optimized path for orthogonalized HTD
        B_root = htd.B{1};  % size: [1 x 1 x r]
        nrm = norm(reshape(B_root, [], 1));  % vector 2-norm
    end
end
