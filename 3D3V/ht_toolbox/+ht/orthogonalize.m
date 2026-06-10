function htd_orth = orthogonalize(htd)
%ORTHOGONALIZE Orthogonalize an HT tensor using QR factorization.
%
%   htd_orth = ht.orthogonalize(htd)
%
%   Input:
%       htd - HTD structure with fields .U, .B, .tree
%   Output:
%       htd_orth - orthogonalized HTD

    tree = htd.tree;

    U = htd.U;
    B = htd.B;
    U_orth = cell(1, tree.N_node);
    B_orth = B;
    R = cell(1, tree.N_node);

    % Step 1: Leaf QR
    for d = 1: tree.orders(1)
        node = tree.dim2ind(d);
        [U_orth{node}, R{node}] = qr(U{node}, "econ");
    end

    % Step 2: Bottom-up orthogonalization
    for i = tree.postorder_nonleaf
        l = tree.children(i, 1);
        r = tree.children(i, 2);

        % Apply (R_r ⊗ R_l) * B
        B_orth{i} = ht.ttm(B_orth{i},R{l},1); %
        B_orth{i} = ht.ttm(B_orth{i},R{r},2); %

        if i ~= 1
            [Q, R{i}] = qr(reshape(B_orth{i},[],htd.rank(i)), "econ");
            B_orth{i} = reshape(Q, size(R{l},1), size(R{r},1), []);
        end
    end

    htd_orth.tree = tree;
    htd_orth.U = U_orth;
    htd_orth.B = B_orth;
    htd_orth.rank = ht.rank(htd_orth);
end
