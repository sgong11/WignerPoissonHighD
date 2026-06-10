function r = rank(x)
%RANK Compute the rank vector from a given HTD structure.
%
%   r = ht.rank(x)
%
%   Input:
%       x - HTD structure with fields .U, .B, .tree, .rank
%
%   Output:
%       r - [1 x N] row vector, where r(i) is the rank at node i

    r = zeros(1, x.tree.N_node);  % Row vector

    % Assign leaf node ranks
    for d = 1: x.tree.orders(1)
        node = x.tree.dim2ind(d);
        r(node) = size(x.U{node}, 2);
    end

    % Assign nonleaf node ranks
    r(x.tree.nonleaf) = cellfun(@(B) size(B, 3), x.B(x.tree.nonleaf));
end