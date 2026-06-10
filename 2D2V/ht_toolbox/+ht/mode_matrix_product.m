function htd_new = mode_matrix_product(htd, A_list, modes)
%MODE_MATRIX_PRODUCT Performs multi-mode matrix product on an HTD.
%
%   htd_new = ht.mode_matrix_product(htd, A_list, modes)
%
%   Input:
%       htd     - HTD structure with fields .U, .B, .tree, .rank
%       A_list  - cell array of matrices {A1, A2, ..., Ak}, each Aμ ∈ ℝ^{m_μ × n_μ}
%       modes   - vector of mode indices μ ∈ [1, d], specifying which modes to transform
%   Output:
%       htd_new - HTD after applying Aμ to Uμ, i.e., Uμ ← Aμ * Uμ for each μ ∈ modes

    htd_new = htd;             % shallow copy to preserve structure
    tree    = htd.tree;

    for i = 1:length(modes)
        mu = modes(i);                 % mode being transformed
        A  = A_list{i};                % transformation matrix Aμ
        node = tree.dim2ind(mu);       % node corresponding to mode-μ

        % Dimension check
        if size(A, 2) ~= size(htd.U{node}, 1)
            error('Dimension mismatch: A{%d} must have %d columns.', ...
                i, size(htd.U{node}, 1));
        end

        % Apply mode-matrix product: Uμ ← Aμ * Uμ
        htd_new.U{node} = A * htd.U{node};

        % Update size and freedom
        htd_new.tree.modesizes(mu) = size(A,1);
        htd_new.tree.freedom(node) = size(A,1);
    end

    % Update internal node freedoms
    for i = tree.postorder_nonleaf
        c = tree.children(i, :);
        htd_new.tree.freedom(i) = htd_new.tree.freedom(c(1)) * htd_new.tree.freedom(c(2));
    end
end
