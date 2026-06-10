% EXAMPLE_GET_SUBTREE_POOL.M
% Demonstrates how to extract all left/right subtrees from an HT tree
% using ht.get_subtree_pool, and how to inspect their mapping relationships.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% Example: a special 6D dimension tree (binary, full, 11 nodes):
%  level 0                   (1)
%                          /     \
%  level 1              (2)       (3)
%                      /  \       /  \
%  level 2          (4)  (5)     (6)  (7)
%                  / \   / \      |    | 
%  level 3       (8)(9) (10)(11)  |    |        
%                 |  |   |   |    |    |
%  modes          1  2   3   4    5    6
children = [2, 3; 
            4, 5; 
            6, 7; 
            8, 9; 
            10, 11; 
            0, 0; 
            0, 0; 
            0, 0; 
            0, 0; 
            0, 0; 
            0, 0];
dim2ind   = [8, 9, 10, 11, 6, 7];
modesizes = [64 64 64 64 64 64];  % each dimension has size 64

% Build tree
tree = ht.make_tree(children, dim2ind, modesizes);
ht.validate_tree(tree);

% Extract all subtrees and mapping pools
[subtree_pool, par2lsub_pool, lsub2par_pool, par2rsub_pool, rsub2par_pool] = ht.get_subtree_pool(tree);

% Display summaries
fprintf('\n========== Subtree Summary ==========\n');
for i = 1:tree.N_node
    if isempty(subtree_pool{i})
        continue;  % skip leaf nodes
    end
    fprintf('Subtree rooted at node %d:\n', i);
    fprintf('  Depth      : %d\n', subtree_pool{i}.depth);
    fprintf('  Num nodes  : %d\n', subtree_pool{i}.N_node);
    fprintf('  Orders     : %s\n', mat2str(subtree_pool{i}.orders));
    
    % Show parent-to-left-subtree map if applicable
    if ~isempty(par2lsub_pool{i})
        nonzero_map = find(par2lsub_pool{i});
        fprintf('  Parent → Left Subtree mapping: %s\n', mat2str(nonzero_map));
    end
    if ~isempty(par2rsub_pool{i})
        nonzero_map = find(par2rsub_pool{i});
        fprintf('  Parent → Right Subtree mapping: %s\n', mat2str(nonzero_map));
    end
    fprintf('\n');
end
