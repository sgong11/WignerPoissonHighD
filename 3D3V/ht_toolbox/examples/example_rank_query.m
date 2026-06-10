% EXAMPLE_RANK_QUERY.M
% Demonstrates how to extract the rank structure from an HT tensor.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% Build tree (same 3D example)
children = [2,3; 4,5; 0,0; 0,0; 0,0];
dim2ind  = [3, 4, 5];
modesizes = [10, 20, 30];
tree = ht.make_tree(children, dim2ind, modesizes);

% Build random HTD
r = [1 3 2 1 1];  % custom rank pattern
A = ht.rand(tree, r);

% Get rank vector
r_query = ht.rank(A);

% Display result
disp('Queried rank vector:');
disp(r_query);
