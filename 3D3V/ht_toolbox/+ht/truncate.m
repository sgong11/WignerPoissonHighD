function htd_new = truncate(htd,tolerance,min_rank,max_rank)
%TRUNCATE Truncate an HTD tensor using Gramian-based SVD thresholding
%
%   htd_new = ht.truncate(htd, tolerance, min_rank, max_rank)
%
%   This function compresses a hierarchical Tucker decomposition (HTD)
%   using a Frobenius norm-based truncation strategy guided by local
%   Gramians. Ranks are selected adaptively within a prescribed range.
%
%   INPUT:
%     htd        - HTD structure with fields .U, .B, .tree, .rank
%     tolerance  - global Frobenius truncation error tolerance
%     min_rank   - vector of minimum allowed ranks at each node
%     max_rank   - vector of maximum allowed ranks at each node
%
%   OUTPUT:
%     htd_new    - compressed HTD tensor with updated ranks

    % Compute local tolerance based on root order
    tol_local = tolerance/sqrt(2*htd.tree.orders(1)-2);

    % Orthogonalize input HTD to prepare for Gramian computation
    htd_new = ht.orthogonalize(htd);

    % Compute Gramians based on orthogonalized HTD structure
    G = ht.gramians_orthog(htd_new);

    % Initialize left singular matrix container of gramians for all nodes
    U_G = cell(1,htd.tree.N_node);

    % === Step 1: Truncate leaf frames ===
    for d = 1: htd.tree.orders(1)
        node = htd.tree.dim2ind(d);

        % SVD on symmetrized Gramian to compute singular values
        [U_G{node}, S] = svd((G{node}+G{node}')/2);
        s_G            = diag(S);
        s_G(find(s_G<s_G(1)*1e-15)) = 0;
        s_U = sqrt(s_G);

        % Choose rank based on Frobenius criterion
        rk  = ht.choose_rank_fro(s_U, tol_local, min_rank(node), max_rank(node));

        U_G{node}          = U_G{node}(:,1:rk);
        htd_new.U{node}    = htd_new.U{node}*U_G{node};
        htd_new.rank(node) = rk;
    end

    % === Step 2: Truncate internal transfer tensors (bottom-up) ===
    for i = htd.tree.postorder_nonleaf(1:end-1)
        % SVD of Gramian to determine optimal rank
        [U_G{i}, S] = svd((G{i}+G{i}')/2);

        s_G = diag(S);
        s_G(find(s_G<s_G(1)*1e-15)) = 0;
        s_U = sqrt(s_G);
        rk  = ht.choose_rank_fro(s_U, tol_local, min_rank(i), max_rank(i));

        % Child indices
        l = htd.tree.children(i,1);
        r = htd.tree.children(i,2);

        % Form new transfer tensors from bottom-up
        U_G{i} = U_G{i}(:,1:rk);
        B = ht.ttm(htd_new.B{i},U_G{i}.',3);
        B = ht.ttm(B,U_G{l}',1);
        htd_new.B{i} = ht.ttm(B,U_G{r}',2);
        htd_new.rank(i) = rk;
    end

    % === Step 3: Form root transor 'matrix' ===
    l = htd.tree.children(1,1);
    r = htd.tree.children(1,2);
    B = ht.ttm(htd_new.B{1},U_G{l}',1);
    htd_new.B{1} = ht.ttm(B,U_G{r}',2);

end