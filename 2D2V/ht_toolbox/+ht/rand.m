function htd = rand(tree, rank)
%RAND Generate a random HTD structure with specified rank.
%
%   htd = ht.rand(tree, rank)
%
%   Inputs:
%       tree - a tree struct from ht.make_tree(...)
%       rank - [1 x N] rank vector (same length as tree.children)
%
%   Output:
%       htd - struct with fields:
%             .tree   : input tree
%             .U      : leaf basis matrices
%             .B      : transfer tensors at internal nodes
%             .rank   : copied from input

    htd.tree = tree;
    htd.rank = rank;
    htd.U = cell(1, tree.N_node);
    htd.B = cell(1, tree.N_node);

    % Assign leaf basis matrices
    for d = 1: tree.orders(1)
        node = tree.dim2ind(d);
        r = rank(node);
        n = tree.modesizes(d);
        htd.U{node} = randn(n, r);  % standard Gaussian
    end

    % Assign transfer tensors at nonleaf nodes
    for i = tree.postorder_nonleaf
        c = tree.children(i, :);
        r_left  = rank(c(1));
        r_right = rank(c(2));
        r_here  = rank(i);
        htd.B{i} = randn(r_left, r_right, r_here);
    end
end