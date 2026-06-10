% EXAMPLE_ADD_AND_SCALAR.M
% Demonstrates HT tensor addition, scalar multiplication, and scalar addition.
% Validates correctness by evaluating at 100 random multi-indices.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% Build tree for a 3D tensor
children = [2 3; 4 5; 0 0; 0 0; 0 0];
dim2ind  = [3 4 5];
modesizes = [10 10 10];
tree = ht.make_tree(children, dim2ind, modesizes);

% Construct more diverse rank pattern: root is 1, others ∈ [1, N]
N = min(modesizes);
rank = ones(1, tree.N_node);
rank(2:end) = randi(N, 1, tree.N_node - 1);

% Generate two random HT tensors
A = ht.rand(tree, rank);
B = ht.rand(tree, rank);

% Combine via addition and scalar operations
C = ht.add(A, B); % A and B must have the same tree structure
A2 = ht.multiply_scalar(A, 2);
A_shifted = ht.add_scalar(A, 5);

% Generate 100 random index sets in each mode
Q = 100;
rng(0);  % for reproducibility
id1 = randi(modesizes(1), Q, 1);
id2 = randi(modesizes(2), Q, 1);
id3 = randi(modesizes(3), Q, 1);

% Evaluate values at 100 indices
vA  = ht.evaluate_index(A, id1, id2, id3);
vB  = ht.evaluate_index(B, id1, id2, id3);
vC  = ht.evaluate_index(C, id1, id2, id3);       % should ≈ vA + vB
vA2 = ht.evaluate_index(A2, id1, id2, id3);      % should ≈ 2 * vA
vAs = ht.evaluate_index(A_shifted, id1, id2, id3);  % should ≈ vA + 5

% Display relative errors
fprintf('--- Value comparison on 100 random points ---\n');
fprintf('‖C - (A + B)‖ / ‖C‖     = %.2e\n', norm(vC - (vA + vB)) / norm(vC));
fprintf('‖A2 - 2A‖ / ‖A2‖        = %.2e\n', norm(vA2 - 2*vA) / norm(vA2));
fprintf('‖As - (A + 5)‖ / ‖As‖   = %.2e\n', norm(vAs - (vA + 5)) / norm(vAs));
