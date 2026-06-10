% EXAMPLE_SQUEEZE.M
% Demonstrates how to use ht.squeeze to remove singleton modes
% and restructure the HT tree accordingly.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% Construct a 6D tree where some physical dimensions are singleton
% Here, we make dims 1, 3, and 6 have mode size 1
modesizes = [1, 32, 1, 64, 48, 1];

% Define a binary full tree for 6D
children = [2, 3;
            4, 5;
            6, 7;
            8, 9;
            10,11;
            0, 0;
            0, 0;
            0, 0;
            0, 0;
            0, 0;
            0, 0];
dim2ind = [8, 9, 10, 11, 6, 7];  % leaf node indices for dimensions 1-6

tree = ht.make_tree(children, dim2ind, modesizes);

% Generate a random HT tensor based on this tree
rank = ones(1, tree.N_node);
rank(2:end) = randi([1, 5], 1, tree.N_node-1);
A = ht.rand(tree, rank);

% Display original tree info
disp('--- Original HTD ---');
disp(['Modesizes : ', mat2str(tree.modesizes)]);
disp(['Tree depth: ', num2str(tree.depth)]);
disp(['Orders    : ', mat2str(tree.orders)]);

% Compute norm before squeeze
norm_before = ht.norm(A);

% Squeeze the HTD to remove singleton modes
A_sq = ht.squeeze(A);

% Display squeezed tree info
disp('--- Squeezed HTD ---');
disp(['Modesizes : ', mat2str(A_sq.tree.modesizes)]);
disp(['Tree depth: ', num2str(A_sq.tree.depth)]);
disp(['Orders    : ', mat2str(A_sq.tree.orders)]);

% Compare norms to confirm consistency
norm_after = ht.norm(A_sq);
rel_err = abs(norm_before - norm_after) / norm_before;

disp(['Original Frobenius norm: ', num2str(norm_before)]);
disp(['Squeezed Frobenius norm: ', num2str(norm_after)]);
disp(['Relative error         : ', num2str(rel_err)]);
