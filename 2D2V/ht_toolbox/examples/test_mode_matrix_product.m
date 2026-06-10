
% test_mode_matrix_product_apply_funciton.m - test ht.mode_matrix_product
% and ht.mode_apply_function
toolbox_root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(toolbox_root);

addpath('htucker_1.2\');

fprintf("===== Test: 4D tensor (balanced tree) =====\n");

% --- 4D HTD with balanced tree ---
N_perdim = 64;
modesizes4 = N_perdim*ones(1,4);
rank4 = [1 20 20 20 30 20 30];

% Balanced binary tree:
children4 = [2, 3; 4, 5; 6, 7; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind4  = [4, 5, 6, 7];  % map mode i → node index
tree4 = ht.make_tree(children4, dim2ind4, modesizes4);
x4 = ht.rand(tree4, rank4);

x4_ref = htensor(tree4.children,tree4.dim2ind,x4.U,x4.B,false);

A2 = randn(10, modesizes4(2));  % new size along mode-2
A4 = randn(5,  modesizes4(4));  % new size along mode-4

tic
for k = 1: 1000
    x_new = ht.mode_matrix_product(x4, {A2 A4}, [2 4]);
end
toc

tic
for k = 1: 1000
    x_ref = ttm(x4_ref, {A2 A4}, [2 4]);
end
toc

x_new_ht = htensor(x_new.tree.children,x_new.tree.dim2ind,x_new.U,x_new.B,false);

fprintf("Reletive diff (4D): %.2e\n", norm(orthog(x_new_ht-x_ref))/norm(x_ref));

fprintf("===== Test: 6D tensor (custom tree) =====\n");

% --- 6D HTD with custom tree ---
% N_perdim = 128;
modesizes6 = N_perdim*ones(1,6);
rank6 = [1 10 10 10 10 20 30 20 30 20 30];
children6 = [2, 3; 4, 5; 6, 7; 8, 9; 10, 11; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind6 = [8, 9, 10, 11, 6, 7];
tree6 = ht.make_tree(children6, dim2ind6, modesizes6);
x6 = ht.rand(tree6, rank6);

x6_ref = htensor(tree6.children,tree6.dim2ind,x6.U,x6.B,false);

A2 = randn(10, modesizes6(2));  % new size along mode-2
A4 = randn(5,  modesizes6(4));  % new size along mode-4
A5 = randn(7,  modesizes6(5));  % new size along mode-4

tic
for k = 1: 1000
    x_new = ht.mode_matrix_product(x6, {A2 A4 A5}, [2 4 5]);
end
toc

tic
for k = 1: 1000
    x_ref = ttm(x6_ref, {A2 A4 A5}, [2 4 5]);
end
toc

x_new_ht = htensor(x_new.tree.children,x_new.tree.dim2ind,x_new.U,x_new.B,false);

fprintf("Reletive diff (4D): %.2e\n", norm(orthog(x_new_ht-x_ref))/norm(x_ref));
