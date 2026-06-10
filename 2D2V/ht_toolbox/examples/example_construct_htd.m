% EXAMPLE_CONSTRUCT_HTD.M
% Demonstrates how to construct HT tensors with all-zero, all-one,
% random real, and random complex entries.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% Build a simple 3D binary HT tree: ((1 ⊗ 2) ⊗ 3)
children = [2,3; 4,5; 0,0; 0,0; 0,0];
dim2ind  = [3, 4, 5];
modesizes = [6, 6, 6];
tree = ht.make_tree(children, dim2ind, modesizes);

% Construct HT tensors
A0 = ht.zeros(tree);                     % all zeros
A1 = ht.ones(tree);                      % all ones
r    = 2 * ones(1, tree.N_node);         % uniform rank-2
r(1) = 1;                                % root rank is fixed to 1   
Ar = ht.rand(tree, r);                   % real random
Ac = ht.rand_complex(tree, r);           % complex random

% Check and display ranks
disp('Ranks of constructed HTDs:');
disp(['zeros     : ', mat2str(ht.rank(A0))]);
disp(['ones      : ', mat2str(ht.rank(A1))]);
disp(['rand      : ', mat2str(ht.rank(Ar))]);
disp(['rand_cplx : ', mat2str(ht.rank(Ac))]);
