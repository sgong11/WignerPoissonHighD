% EXAMPLE_MODE_TRANSFORM.M
% Demonstrates how to transform specific modes of an HT tensor
% using mode-matrix product and elementwise function application.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% Build a 3D HT tree
children = [2 3; 4 5; 0 0; 0 0; 0 0];
dim2ind  = [3 4 5];
modesizes = [12 10 8];  % different sizes to test dim changes
tree = ht.make_tree(children, dim2ind, modesizes);

% Random rank with root = 1
N = min(modesizes);
rank = ones(1, tree.N_node);
rank(2:end) = randi(N, 1, tree.N_node - 1);

% Create random HT tensor
A = ht.rand(tree, rank);

% === Mode-matrix product ===
% Reduce mode-1 (12→6) by applying a 6x12 matrix to U₁
A1 = ht.mode_matrix_product(A, {randn(6,12)}, [1]);
fprintf('Mode-matrix product: mode-1 size changed from 12 to %d\n', ...
    A1.tree.modesizes(1));

% === Mode-function application ===
% Apply elementwise nonlinear function to mode-2 (e.g., abs)
F_list = {@(x) abs(x)};
A2 = ht.mode_apply_function(A, F_list, [2]);

% Display size change and quick norm check
fprintf('Function application: norm(A) = %.4e, norm(A2) = %.4e\n', ...
    ht.norm(A), ht.norm(A2));

% === Multiple-mode transformation (matrix product and function) ===

% Mode 1: 12 → 6 via matrix
% Mode 3: 8 → 4 via matrix
A_multi = ht.mode_matrix_product(A, {randn(6,12), randn(4,8)}, [1, 3]);
fprintf('\n--- Multiple-mode Matrix Product ---\n');
fprintf('New mode sizes: [%d, %d, %d]\n', A_multi.tree.modesizes);

% Mode 2 & 3: apply nonlinear function (e.g., square)
F2 = {@(x) x.^2, @(x) sin(x)};
A_func_multi = ht.mode_apply_function(A, F2, [2, 3]);
fprintf('\n--- Multiple-mode Function Application ---\n');
fprintf('norm(A) = %.4e, norm(A_func_multi) = %.4e\n', ...
    ht.norm(A), ht.norm(A_func_multi));
