% EXAMPLE_NORM_AND_INNER.M
% Demonstrates Frobenius norm and inner product of HT tensors.
% Shows the importance of orthogonalization for stable norm computation.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% Build 3D tree
children = [2 3; 4 5; 0 0; 0 0; 0 0];
dim2ind  = [3 4 5];
modesizes = [10 10 10];
tree = ht.make_tree(children, dim2ind, modesizes);

% Random rank (root=1, others random)
N = min(modesizes);
rank = ones(1, tree.N_node);
rank(2:end) = randi(N, 1, tree.N_node - 1);

% Construct random HT tensor
A = ht.rand(tree, rank);

% Compute norm and inner product
nrm_A  = ht.norm(A);
ip_AA  = ht.inner_product(A, A);

% Display result
fprintf('--- Norm and Inner Product ---\n');
fprintf('‖A‖_F                 = %.6e\n', nrm_A);
fprintf('⟨A, A⟩                = %.6e ≈ ‖A‖² = %.6e\n', ip_AA, nrm_A^2);


% NOTE:
% For HT tensors that are not orthogonalized, the Frobenius norm may
% exhibit significant numerical instability due to the recursive contraction
% structure. For example, evaluating ‖A - A‖ without orthogonalization may
% yield a relative norm on the order of 1e-8 or worse, even though the result
% is theoretically zero.
%
% To ensure accurate norm-based comparisons (e.g., checking ‖A - B‖ ≈ 0),
% always use:
%     A_orth = ht.orthogonalize(A);
%     norm(A_orth, 'orth');

% Check error introduced by subtraction
C = ht.add(A, ht.multiply_scalar(A, -1));  % C ≈ 0 (in theory)
rel_norm_C = ht.norm(C) / nrm_A;           % may be large due to instability

% Orthogonalize before norm computation
C_orth = ht.orthogonalize(C);
rel_norm_C_orth = ht.norm(C_orth, 'orth') / nrm_A;

fprintf('\n--- Numerical Stability Test ---\n');
fprintf('‖A - A‖ / ‖A‖ (raw)         = %.2e\n', rel_norm_C); 
fprintf('‖A - A‖ / ‖A‖ (orthogonal)  = %.2e\n', rel_norm_C_orth);
% if you see ‖A - A‖ / ‖A‖ (raw) = 0.00e+00, run this test for a few more 
% times 