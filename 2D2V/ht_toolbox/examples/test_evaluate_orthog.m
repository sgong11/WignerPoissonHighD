
% test_evaluate.m - test ht.evaluate_index vs ht.evaluate_index_fast
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
x4 = ht.rand(tree4, rank4);

% Generate random batch of indices
Ntest = N_perdim;
id4 = randi([1 64], Ntest, 4);
id4 = mat2cell(id4,Ntest,[1 1 1 1]);

tic
for k = 1: 1000
    val1 = ht.evaluate_index(x4, id4{:});
end
toc

tic
for k = 1: 1000
    val2 = get_val_ht_2D2V(x4,id4{:},x4.rank);
end
toc
fprintf("Max abs diff (4D): %.2e\n", max(abs(val1 - val2))/max(abs(val1)));

tic
for k = 1: 1000
    val1 = ht.evaluate_fiber(x4, 1,1:N_perdim,5,3);
end
toc

% tic
% for k = 1: 1000
%     val2 = ht.evaluate_fiber_new(x4, 1,1:N_perdim,5,3);
% end
% toc

id4_ = [1*ones(N_perdim,1)  (1:N_perdim)' 5*ones(N_perdim,1)  3*ones(N_perdim,1)];
id4_ = mat2cell(id4_,Ntest,[1 1 1 1]);

tic
for k = 1: 1000
    val3 = ht.evaluate_index(x4, id4_{:});
end
toc
fprintf("Max abs diff (4D): %.2e\n", max(abs(val1 - val3))/max(abs(val1)));

x4_orthog = ht.orthogonalize(x4);

for i = 4: 7
    U = x4_orthog.U{i};
    orth_err = norm(U' * U - eye(size(U, 2)));
    fprintf("Node %d orthogonality error: %.2e\n", i, orth_err);
end

val1 = ht.evaluate_index(x4, id4{:});
val2 = ht.evaluate_index(x4_orthog, id4{:});
fprintf("Max abs diff (4D): %.2e\n", max(abs(val1 - val2))/max(abs(val1)));

fprintf("===== Test: 6D tensor (custom tree) =====\n");

% --- 6D HTD with custom tree ---
% N_perdim = 128;
modesizes6 = N_perdim*ones(1,6);
rank6 = [1 10 10 10 10 20 30 20 30 20 30];
children6 = [2, 3; 4, 5; 6, 7; 8, 9; 10, 11; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind6 = [8, 9, 10, 11, 6, 7];
tree6 = ht.make_tree(children6, dim2ind6, modesizes6);
x6 = ht.rand(tree6, rank6);

id6 = randi([1 64], Ntest, 6);
id6 = mat2cell(id6,64,[1 1 1 1 1 1]);

tic
for k = 1: 1000
val1 = ht.evaluate_index(x6, id6{:});
end
toc

tic
for k = 1: 1000
val2 = get_val_ht_3D3V(x6,id6{:},x6.rank);
end
toc

fprintf("Max abs diff (6D): %.2e\n", max(abs(val1 - val2))/max(abs(val1)));

x6_orthog = ht.orthogonalize(x6);

for i = 1: 6
    id = x6.tree.dim2ind(i);
    U = x6_orthog.U{id};
    orth_err = norm(U' * U - eye(size(U, 2)));
    fprintf("Node %d orthogonality error: %.2e\n", i, orth_err);
end

val1 = ht.evaluate_index(x6, id6{:});
val2 = ht.evaluate_index(x6_orthog, id6{:});
fprintf("Max abs diff (6D): %.2e\n", max(abs(val1 - val2))/max(abs(val1)));

tic
for k = 1: 1000
    val1 = ht.evaluate_fiber(x6, 1,1:N_perdim,5,3,20,30);
end
toc

% tic
% for k = 1: 1000
%     val2 = ht.evaluate_fiber_new(x6, 1,1:N_perdim,5,3,20,30);
% end
% toc

id6_ = [1*ones(N_perdim,1)  (1:N_perdim)' 5*ones(N_perdim,1)  3*ones(N_perdim,1) 20*ones(N_perdim,1) 30*ones(N_perdim,1)];
id6_ = mat2cell(id6_,Ntest,[1 1 1 1 1 1]);

tic
for k = 1: 1000
    val3 = ht.evaluate_index(x6, id6_{:});
end
toc

tic
for k = 1: 1000
    val4 = get_fiber1or2_ht_3D3V(x6, 1,1:N_perdim,5,3,20,30, x6.rank, Ntest);
end
toc
fprintf("Max abs diff (6D): %.2e\n", max(abs(val1 - val3))/max(abs(val1)));

