% EXAMPLE_BUILD_TREE.M
% Demonstrates how to construct a Hierarchical Tucker dimension tree
% using ht.make_tree, validate its structure, and extract subtrees.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

% Tree is defined by the 'children' array:
%   - children(i,:) = [l, r] means node i has left child l and right child r.
%   - children(i,:) = [0, 0] means node i is a leaf node (i.e., corresponds to a physical dimension).
%
% The physical dimensions are mapped to leaf node indices by 'dim2ind'
%
% modesizes is not a relevant property of tree. We encode this property for
% convenience of the ht.HTACA function

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

% Build the tree using ht package
tree = ht.make_tree(children, dim2ind, modesizes);

% Validate the tree
ht.validate_tree(tree);

% Display tree structure
disp('--- Tree Summary ---');
disp(['Number of nodes      : ', num2str(tree.N_node)]);
disp(['Tree depth           : ', num2str(tree.depth)]);
disp(['Postorder nonleaf    : ', mat2str(tree.postorder_nonleaf)]);
disp(['Orders               : ', mat2str(tree.orders)]);
disp(['Freedom (mode sizes) : ', mat2str(tree.freedom)]);

% Visualize subtree extraction
[left_subtree, ~, ~] = ht.get_left_subtree(tree);
[right_subtree, ~, ~] = ht.get_right_subtree(tree);

disp('--- Left Subtree Summary ---');
disp(['Depth: ', num2str(left_subtree.depth), ...
      ', Nodes: ', num2str(left_subtree.N_node), ...
      ', Orders: ', mat2str(left_subtree.orders)]);

disp('--- Right Subtree Summary ---');
disp(['Depth: ', num2str(right_subtree.depth), ...
      ', Nodes: ', num2str(right_subtree.N_node), ...
      ', Orders: ', mat2str(right_subtree.orders)]);
