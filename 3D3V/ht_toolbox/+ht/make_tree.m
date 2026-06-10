function tree = make_tree(children, dim2ind, modesizes)
%MAKE_TREE Constructs a standardized tree structure.
%
%   tree = ht.make_tree(children, dim2ind, modesizes)
%
%   Inputs:
%       children   - [N x 2] matrix, each row is [left_child, right_child]
%       dim2ind    - [1 x D] vector, dim2ind(d) = index of node representing mode d
%       modesizes  - [1 x D] vector, size of each mode
%
%   Output:
%       tree - struct with fields:
%         .children               : [N x 2] matrix
%         .dim2ind                : [1 x D] vector
%         .modesizes              : [1 x D] vector
%         .nonleaf                : [1 x N] logical (true if node is nonleaf)
%         .tree.postorder_nonleaf : [1 x d-1] post-order traversal for nonlieaf nodes
%         .orders                 : [1 x N] number of modes per node
%         .freedom                : [1 x N] product of mode sizes per node
%         .depth                  : scalar, max tree depth (edge count)
%         .N_node                 : scalar, number of nodes (including root)

    % Validate dimensions
    N = size(children, 1);
    D = length(modesizes);

    if length(dim2ind) ~= D
        error('Length of dim2ind and modesizes must match.');
    end

    % Initialize tree
    tree.children  = children;
    tree.dim2ind   = dim2ind(:)';  % row vector
    tree.modesizes = modesizes(:)';

    % Compute derived fields
    tree.nonleaf = (children(:,1) ~= 0) | (children(:,2) ~= 0);  % non-leaf nodes
    tree.nonleaf = tree.nonleaf(:)';

    % Number of nodes
    tree.N_node = N;
    
    % Computes post-order traversal for nonleaf nodes
    tree.postorder_nonleaf = ht.postorder_nonleaf(tree);

    % Initialize orders/freedom
    tree.orders  = zeros(1, N);
    tree.freedom = ones(1, N);

    for d = 1: D
        node = dim2ind(d);
        tree.orders(node)  = 1;
        tree.freedom(node) = modesizes(d);
    end

    % Post-order traversal for orders and freedom
    for i = tree.postorder_nonleaf
        c = children(i, :);
        tree.orders(i)  = tree.orders(c(1)) + tree.orders(c(2));
        tree.freedom(i) = tree.freedom(c(1)) * tree.freedom(c(2));
    end

    % Compute depth
    tree.depth = ht.tree_depth(tree.children);


end
