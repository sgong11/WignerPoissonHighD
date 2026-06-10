function min_rank = get_default_min_rank(tree)
%GET_DEFAULT_MIN_RANK Return the default minimum rank for each node in the tree
%
%   min_rank = ht.get_default_min_rank(tree)
%
%   This function initializes a minimal rank of 1 for every node
%   in a hierarchical Tucker tree. It serves as a conservative lower
%   bound for rank truncation during HTD compression.

    % Set all node ranks (leaves and nonleaves) to a minimum of 1
    min_rank = ones(1,tree.N_node);

end