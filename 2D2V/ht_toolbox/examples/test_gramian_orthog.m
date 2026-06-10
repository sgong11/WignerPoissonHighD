% test_gramian_orthog.m - test ht.gramian_orthog
toolbox_root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(toolbox_root);
addpath('htucker_1.2')
clear

fprintf("===== Test: 4D tensor (balanced tree) =====\n");

% --- 4D HTD with balanced tree ---
N_perdim = 64;
modesizes4 = N_perdim*ones(1,4);
rank4 = [1 20 20 20 30 20 30];

% Balanced binary tree:
children4 = [2, 3; 4, 5; 6, 7; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind4  = [4, 5, 6, 7];  % map mode i → node index
tree4 = ht.make_tree(children4, dim2ind4, modesizes4);
x4 = ht.rand_complex(tree4, rank4);
x4_2 = ht.rand_complex(tree4, rank4);


x4_htucker   = htensor(children4,dim2ind4,x4.U,x4.B,false);
x4_2_htucker = htensor(children4,dim2ind4,x4_2.U,x4_2.B,false);

opts.abs_eps  = 1e-6
opts.max_rank = 64;
x_cell = {x4_htucker,x4_2_htucker};
x  = htensor.truncate_sum(x_cell,opts)

x4 = ht.orthogonalize(x4);

x4_htucker = orthog(x4_htucker);

G_1 = ht.gramians_orthog(x4);

G_2 = gramians(x4_htucker);

N_node = size(children4,1);
for i = 1: N_node
    norm(G_1{i} - G_2{i},'fro')/norm(G_1{i},'fro')
end