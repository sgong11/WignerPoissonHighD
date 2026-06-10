function d = tree_depth(children)
%TREE_DEPTH Computes the maximum depth of a subtree.
%
%   d = ht.tree_depth(children)
%
%   Input:
%       children - [n x 2] matrix of local subtree connectivity
%
%   Output:
%       d        - scalar integer, depth (number of edges from root to leaf)

    function depth = dfs(node)
        if all(children(node, :) == 0)
            depth = 0;  % leaf node has depth 0
        else
            depth = 1 + max(dfs(children(node, 1)), dfs(children(node, 2)));
        end
    end

    d = dfs(1);  % root is node 1 in all subtree structures
end