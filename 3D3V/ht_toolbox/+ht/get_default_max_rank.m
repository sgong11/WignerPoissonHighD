function max_rank = get_default_max_rank(tree)
%GET_DEFAULT_MAX_RANK Return the default maximum rank for each node in the tree
%
%   max_rank = ht.get_default_max_rank(tree)
%
%   This function estimates an upper bound for each node's rank, such that it matches the full rank of the node's matricization.
%
%   OUTPUT:
%     max_rank - 1×N vector, where N is the number of tree nodes. Each
%                entry gives the maximum rank allowed at that node.

    % Initialize all ranks to 1 (will be updated below)
    max_rank = ones(1,tree.N_node); % max rank for root is always one
    
    % Set upper bound for leaf nodes
    for d = 1: tree.orders(1)
        % Upper bound is the minimum of:
        % - local freedom
        % - ratio of root freedom to local freedom (for balance)
        node = tree.dim2ind(d);
        max_rank(node) = min(tree.freedom(node),tree.freedom(1)/tree.freedom(node));
    end
    
    % Set upper bound for internal (non-root) nodes (postorder traversal)
    for i = tree.postorder_nonleaf(1:end-1)
        % Same logic
        max_rank(i) = min(tree.freedom(i),tree.freedom(1)/tree.freedom(i));
    end

end