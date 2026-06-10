function order = postorder_nonleaf(tree)
%POSTORDER_NONLEAF Computes post-order traversal for non-leaf nodes.
%
%   order = ht.postorder_nonleaf(tree)
%
%   Input:
%       tree - a struct with fields:
%              .children : [N x 2] array of child indices (0 for leaf)
%              .nonleaf  : logical vector of length N indicating non-leaf nodes
%
%   Output:
%       order - row vector of indices of non-leaf nodes in post-order,
%               i.e., each node appears after its children.

    visited = false(1, tree.N_node);
    order = [];

    function dfs(i)
        if i == 0 || visited(i)
            return;
        end
        dfs(tree.children(i, 2));
        dfs(tree.children(i, 1));
        if tree.nonleaf(i)
            order(end+1) = i; %#ok<AGROW>
        end
        visited(i) = true;
    end

    dfs(1);  % root is always node 1
end
