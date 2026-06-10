% EXAMPLE_EVALUATE_SAMPLES.M
% Demonstrates how to evaluate values from an HT tensor using:
% - evaluate_index
% - evaluate_index_sum
% - evaluate_fiber
% - evaluate_fiber_sum
% - evaluate_slice

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

%% Build a 4D HT tensor and a list of scaled copies

% Dimension tree: binary, 4 leaf nodes (4D tensor)
children = [2 3; 4 5; 6 7; 0 0; 0 0; 0 0; 0 0];
dim2ind = [4 5 6 7];
modesizes = [12 10 8 6];
tree = ht.make_tree(children, dim2ind, modesizes);

% Random rank vector (root = 1, others = random integers)
rank_vec = ones(1, tree.N_node);
rank_vec(2:end) = randi(4, 1, tree.N_node - 1);

% Generate base HT tensor and scaled copies
A = ht.rand(tree, rank_vec);
A_list = {
    A, ...
    ht.multiply_scalar(A, 2), ...
    ht.multiply_scalar(A, -0.5)
};

%% Evaluate multiple scalar entries using evaluate_index

Q = 100;
id1 = randi(modesizes(1), Q, 1);
id2 = randi(modesizes(2), Q, 1);
id3 = randi(modesizes(3), Q, 1);
id4 = randi(modesizes(4), Q, 1);

valA = ht.evaluate_index(A, id1, id2, id3, id4);
fprintf('evaluate_index: evaluated %d random entries.\n', Q);

%% Evaluate sum across multiple tensors using evaluate_index_sum

valSum = ht.evaluate_index_sum(A_list, id1, id2, id3, id4);
fprintf('evaluate_index_sum: max(|sum - expected|) = %.2e\n', ...
    max(abs(valSum - 2.5 * valA)));  % since total = A * (1 + 2 - 0.5) = 2.5A

%% Extract a fiber along mode-2 using evaluate_fiber

% Fix modes 1, 3, 4; vary mode-2
idx = {3, 1:10, 2, 4};
fiber_val = ht.evaluate_fiber(A, idx{:});
fprintf('evaluate_fiber: extracted fiber of length %d (along mode-2)\n', ...
    length(fiber_val));

%% Evaluate sum of fibers across A_list using evaluate_fiber_sum

fiber_val_sum = ht.evaluate_fiber_sum(A_list, idx{:});
fprintf('evaluate_fiber_sum: max(|sum - expected|) = %.2e\n', ...
    max(abs(fiber_val_sum - 2.5 * fiber_val)));

%% Extract a 2D slice along mode-2 and mode-3 using evaluate_slice

% Fix mode-1 and mode-4; vary mode-2 and mode-3
slice_idx = {2, 1:5, 1:6, 4};
slice_val = ht.evaluate_slice(A, slice_idx{:});
fprintf('evaluate_slice: slice size = %d x %d (modes 2 and 3)\n', ...
    size(slice_val, 1), size(slice_val, 2));

%% Notes:
% - All tensors used in this example share the same dimension tree.
% - *_sum variants perform evaluations across multiple HT tensors.
% - evaluate_index: returns scalar values at specified index tuples.
% - evaluate_fiber: returns a 1D array by fixing all but one mode.
% - evaluate_slice: returns a 2D matrix by fixing all but two modes.
