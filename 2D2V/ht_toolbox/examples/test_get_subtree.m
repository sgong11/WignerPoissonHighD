toolbox_root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(toolbox_root);


% Define tree structure
modesizes6 = [64 64 64 64 64 64];
children6 = [2, 3; 4, 5; 6, 7; 8, 9; 10, 11; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind6 = [8, 9, 10, 11, 6, 7];
tree = ht.make_tree(children6, dim2ind6, modesizes6);

[subtree_pool, par2lsub_pool, lsub2par_pool, par2rsub_pool, rsub2par_pool] = ht.get_subtree_pool(tree);

% [tree_l, old2new_l, new2old_l] = ht.get_left_subtree(tree)
% [tree_r, old2new_r, new2old_r] = ht.get_right_subtree(tree)
% tree = ht.make_tree([2 3; 0 0; 0 0],[2 3], [64 64])
% [tree_l, old2new_l, new2old_l] = ht.get_left_subtree(tree)
% [tree_r, old2new_r, new2old_r] = ht.get_right_subtree(tree)

% Define tree structure
children = [2 3; 4 5; 6 7; 8 9; 0 0; 10 11; 0 0; 0 0; 0 0; 0 0; 0 0];
dim2ind  = [8, 9, 5, 10, 11, 7];
modesizes = [64 64 64 64 64 64];
tree = ht.make_tree(children, dim2ind, modesizes);

[subtree_pool, par2lsub_pool, lsub2par_pool, par2rsub_pool, rsub2par_pool] = ht.get_subtree_pool(tree);