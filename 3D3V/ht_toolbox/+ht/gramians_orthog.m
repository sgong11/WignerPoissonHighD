function G = gramians_orthog(htd)
%GRAMIANS_ORTHOG Compute the reduced gramians of an orthogonalized HTD
%
%   G = ht.gramians_orthog(htd)
%
%   INPUT:
%     htd - an orthogonalized HTD structure with fields .U, .B, .tree, .rank
%
%   OUTPUT:
%     G - the reduced gramians of htd

% Initialization
G      = cell(1,htd.tree.N_node);
G{1}   = 1; % The reduced gramian of the root node is always one

% Run DFS traversal from root 
dfs(1)

    function dfs(node)
        if htd.tree.nonleaf(node)
            node_l = htd.tree.children(node,1);
            node_r = htd.tree.children(node,2);

            % The reduced Gramians of the children can be computed from the parent node.
            B_temp = ht.ttm(conj(htd.B{node}),G{node},3);
            G{node_r} = reshape(permute(htd.B{node},[1 3 2]),[],htd.rank(node_r)).'*reshape(permute(B_temp,[1 3 2]),[],htd.rank(node_r));
            G{node_l} = reshape(permute(htd.B{node},[2 3 1]),[],htd.rank(node_l)).'*reshape(permute(B_temp,[2 3 1]),[],htd.rank(node_l));

            % Visit child nodes
            dfs(node_l);
            dfs(node_r);
        end
    end

end