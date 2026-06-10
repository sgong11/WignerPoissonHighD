% Test compress_htd_sum.m
toolbox_root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(toolbox_root);

clear

% Define a 4D balanced tree
children4 = [2, 3; 4, 5; 6, 7; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind4  = [4, 5, 6, 7];  % map mode i → node index
modesizes4 = [64 64 64 64];
tree = ht.make_tree(children4, dim2ind4, modesizes4);


% Generate 50 random htds
for k = 1: 20
    rank{k} = ones(1,7);
    rank{k}(2:7) = randi([1 64],1,6);
    x{k} = ht.rand(tree,rank{k});
end

% Add them
tic
z = ht.compress_htd_sum(x,1e-14);
toc

id = randi([1 64], 100, 4);

tic
for k = 1:10
    val1 = ht.evaluate_index(z,id);
end
toc
tic
for k = 1:10
    val2 = ht.evaluate_index_sum(x,id);
end
toc

max(abs(val1-val2))/max(abs(val1))

% --- 6D HTD with custom tree ---
% N_perdim = 128;
modesizes6 = 32*ones(1,6);
children6 = [2, 3; 4, 5; 6, 7; 8, 9; 10, 11; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind6 = [8, 9, 10, 11, 6, 7];
tree = ht.make_tree(children6, dim2ind6, modesizes6);


% Generate 20 random htds
for k = 1: 20
    rank{k} = ones(1,11);
    rank{k}(2:11) = randi([1 32],1,10);
    x{k} = ht.rand(tree,rank{k});
end

% Add them
tic
z = ht.compress_htd_sum(x,1e-14);
toc

id = randi([1 32], 100, 6);

tic
for k = 1:10
    val1 = ht.evaluate_index(z,id);
end
toc
tic
for k = 1:10
    val2 = ht.evaluate_index_sum(x,id);
end
toc

max(abs(val1-val2))/max(abs(val1))

% --- 6D HTD with custom tree ---
% N_perdim = 128;
children = [2 3; 4 5; 6 7; 8 9; 0 0; 10 11; 0 0; 0 0; 0 0; 0 0; 0 0];
dim2ind  = [8, 9, 5, 10, 11, 7];
modesizes = [32 32 32 32 32 32];
tree = ht.make_tree(children, dim2ind, modesizes);


% Generate 20 random htds
for k = 1: 20
    rank{k} = ones(1,11);
    rank{k}(2:11) = randi([1 32],1,10);
    x{k} = ht.rand(tree,rank{k});
end

% Add them
tic
z = ht.compress_htd_sum(x,1e-14);
toc

id = randi([1 32], 100, 6);

tic
for k = 1:10
    val1 = ht.evaluate_index(z,id);
end
toc
tic
for k = 1:10
    val2 = ht.evaluate_index_sum(x,id);
end
toc

max(abs(val1-val2))/max(abs(val1))

