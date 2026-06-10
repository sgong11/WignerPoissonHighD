toolbox_root = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(toolbox_root);


% Define tree structure
modesizes6 = [64 64 64 64 1 1];
children6 = [2, 3; 4, 5; 6, 7; 8, 9; 10, 11; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind6 = [8, 9, 10, 11, 6, 7];
tree = ht.make_tree(children6, dim2ind6, modesizes6);
rank = [1 10 10 10 10 20 30 20 30 20 30];

htd = ht.rand(tree, rank);

htd_new = ht.squeeze(htd);

id4 = randi([1 64], 64, 4);
id6 = [id4(:,1:4) ones(64,1) ones(64,1)];
id4 = mat2cell(id4, 64, [1 1 1 1]);
id6 = mat2cell(id6, 64, [1 1 1 1 1 1]);

val1 = ht.evaluate_index(htd,id6{:});
val2 = ht.evaluate_index(htd_new,id4{:});

max(abs(val1 - val2))/max(val1)

% Define tree structure
modesizes6 = [1 1 1 1 64 64];
children6 = [2, 3; 4, 5; 6, 7; 8, 9; 10, 11; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0];
dim2ind6 = [8, 9, 10, 11, 6, 7];
tree = ht.make_tree(children6, dim2ind6, modesizes6);
rank = [1 10 10 10 10 20 30 20 30 20 30];

htd = ht.rand(tree, rank);

htd_new = ht.squeeze(htd);

id2 = randi([1 64], 64, 2);
id6 = [ones(64,4) id2(:,1:2)];
id2 = mat2cell(id2, 64, [ 1 1]);
id6 = mat2cell(id6, 64, [1 1 1 1 1 1]);

val1 = ht.evaluate_index(htd,id6{:});
val2 = ht.evaluate_index(htd_new,id2{:});

max(abs(val1 - val2))/max(val1)

% Define tree structure
children = [2 3; 4 5; 6 7; 8 9; 0 0; 10 11; 0 0; 0 0; 0 0; 0 0; 0 0];
dim2ind  = [8, 9, 5, 10, 11, 7];
modesizes = [64 64 1 1 64 1];
tree = ht.make_tree(children, dim2ind, modesizes);

htd = ht.rand(tree, rank);

htd_new = ht.squeeze(htd);

id3 = randi([1 64], 64, 3);
id6 = [id3(:,1:2) ones(64,1) ones(64,1) id3(:,3) ones(64,1)];
id3 = mat2cell(id3, 64, [1 1 1]);
id6 = mat2cell(id6, 64, [1 1 1 1 1 1]);

val1 = ht.evaluate_index(htd,id6{:});
val2 = ht.evaluate_index(htd_new,id3{:});

max(abs(val1 - val2))/max(val1)

% Define tree structure
children = [2 3; 4 5; 6 7; 8 9; 0 0; 10 11; 0 0; 0 0; 0 0; 0 0; 0 0];
dim2ind  = [8, 9, 5, 10, 11, 7];
modesizes = [64 1 64 1 64 1];
tree = ht.make_tree(children, dim2ind, modesizes);

htd = ht.rand(tree, rank);

htd_new = ht.squeeze(htd);

id3 = randi([1 64], 100, 3);
id6 = [id3(:,1) ones(100,1) id3(:,2) ones(100,1) id3(:,3) ones(100,1)];
id3 = mat2cell(id3, 100, [1 1 1]);
id6 = mat2cell(id6, 100, [1 1 1 1 1 1]);

val1 = ht.evaluate_index(htd,id6{:});
val2 = ht.evaluate_index(htd_new,id3{:});

max(abs(val1 - val2))/max(val1)