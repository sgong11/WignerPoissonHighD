function htd_new = mode_apply_function(htd, F_list, modes)
%MODE_APPLY_FUNCTION Applies functions to specified modes of HTD leaf bases.
%
%   htd_new = ht.mode_apply_function(htd, F_list, modes)
%
%   Input:
%       htd     - HTD structure with fields .U, .B, .tree, .rank
%       F_list  - cell array of function handles {@f1, @f2, ..., @fk}, each fμ: ℝ^{nμ×rμ} → ℝ^{mμ×rμ}
%       modes   - vector of mode indices μ ∈ [1, d], specifying which modes to transform
%
%   Output:
%       htd_new - HTD after applying fμ to Uμ, i.e., Uμ ← fμ(Uμ) for each μ ∈ modes

    htd_new = htd;  % shallow copy
    tree = htd.tree;

    if isa(F_list, 'function_handle')
        F_list = repmat({F_list}, 1, length(modes));
    end

    for i = 1: length(modes)
        mu = modes(i);
        fmu = F_list{i};
        node = tree.dim2ind(mu);

        U_new = fmu(htd.U{node});  % Apply function handle
        htd_new.U{node} = U_new;

        % Update mode size and freedom
        htd_new.tree.modesizes(mu) = size(U_new, 1);
        htd_new.tree.freedom(node) = size(U_new, 1);
    end

    % Update internal node freedom
    for i = tree.postorder_nonleaf
        c = tree.children(i, :);
        htd_new.tree.freedom(i) = htd_new.tree.freedom(c(1)) * htd_new.tree.freedom(c(2));
    end
end
