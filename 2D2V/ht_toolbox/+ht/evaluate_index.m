function val = evaluate_index(htd, varargin)
%EVALUATE_INDEX Evaluates the HT tensor at given multi-indices.
%
%   val = ht.evaluate_index(htd, id1, id2, ..., idD)
%
%   Input:
%       htd  - HTD structure with fields .U, .B, .tree, .rank
%       id1,...,idD - index arguments
%       each of the indix arguments is a [Q x 1] vector containing Q
%       indices of the corresponding mode
%
%   Output:
%       val  - [Q x 1] vector of values at each index

    if length(varargin) == 1
        % htd is vector
        val = htd.U{varargin{1}};
        return
    end

    tree = htd.tree;
    D   = tree.orders(1);   % number of modes
    Q   = length(varargin{1});           % number of index queries

    % Step 1: collect leaf evaluations
    V = cell(1, tree.N_node);  % intermediate results
    for d = 1: D
        node = tree.dim2ind(d);
        U = htd.U{node};  % [N_d x r_d]
        V{node} = U(varargin{d}, :);  % [Q x r_d]
    end

    % Step 2: evaluate internal nodes bottom-up
    for i = tree.postorder_nonleaf
        l = tree.children(i, 1);
        r = tree.children(i, 2);
        B = htd.B{i};  % [r_l x r_r x r_i]

        % Contract over left and right
        A = V{l};  % [Q x r_l]
        C = V{r};  % [Q x r_r]
        [Q1, r_l] = size(A);
        [Q2, r_r] = size(C);
        r_i = size(B, 3);

        if Q1 ~= Q2 || size(B,1)~=r_l || size(B,2)~=r_r
            error("Dimension mismatch at node %d", i);
        end

        % Evaluate: v_i(q, :) = sum_{a,b} A(q,a) * B(a,b,:) * C(q,b)
        if Q*r_r*r_i < 1e5
            % Use vectorized fast path for small size, fallback to loop if too large
            % Threshold ~1e5 elements ≈ 0.8 MB for double precision
            AB   = pagemtimes(A,reshape(B,r_l,[]));
            CC   = repmat(C, 1, r_i);
            V_   = reshape(AB.*CC,Q,r_r,r_i);
            V{i} = reshape(sum(V_,2),Q,r_i);
        else
            V{i} = zeros(Q, r_i);
            for k = 1:r_i
                Bk = B(:, :, k);
                V{i}(:, k) = sum(A .* (C * Bk.'), 2);
            end
        end
    end

    % Step 3: output from root node (node 1)
    val = V{1};  % [Q x 1]
end
