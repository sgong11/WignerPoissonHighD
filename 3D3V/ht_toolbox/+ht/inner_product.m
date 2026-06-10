function ip = inner_product(htd1, htd2)
%INNER_PRODUCT Computes ⟨X, Y⟩ where X and Y are HT tensors.
%
%   ip = ht.inner_product(htd1, htd2)
%
%   Input:
%       htd1, htd2 - HTD structures with fields .U, .B, .tree, .rank
%
%   Output:
%       ip - scalar value of the inner product ⟨htd1, htd2⟩

    tree = htd1.tree;
    M = cell(1, tree.N_node);  % intermediate contractions

    % Leaf node: Mt = Ux' * Uy
    for d = 1: tree.orders(1)
        node = tree.dim2ind(d);
        Ux = htd1.U{node};   % [n_d x r_d]
        Uy = htd2.U{node};   % [n_d x r_d]
        M{node} = Ux' * Uy;  % [r_d x r_d]
    end

    % Nonleaf node: Mt = Bx^H * (Ml ⊗ Mr) * By
    for i = tree.postorder_nonleaf
        l = tree.children(i, 1);
        r = tree.children(i, 2);
        Blx = htd1.B{i};  % [r_l_1 x r_r_1 x r_i_1]
        Bly = htd2.B{i};  % [r_l_2 x r_r_2 x r_i_2]

        Ml = M{l};  % [r_l_1 x r_l_2]
        Mr = M{r};  % [r_r_1 x r_r_2]

        r_i_1 = htd1.rank(i);
        r_i_2 = htd2.rank(i);

        Mtmp = zeros(r_i_1, r_i_2);

        for a = 1: r_i_1
            Ba = Blx(:, :, a);  % [r_l_1 x r_r_1]
            for b = 1:r_i_2
                Bb = Bly(:, :, b);  % [r_l_2 x r_r_2]
                % Efficient contraction without forming kron(Ml, Mr)
                tmp = Ml * Bb * Mr';
                Mtmp(a, b) = Ba(:)' * tmp(:);
            end
        end

        M{i} = Mtmp;
    end

    ip = M{1};  % ⟨htd1, htd2⟩ ∈ ℝ
end
