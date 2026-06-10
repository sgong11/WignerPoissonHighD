function val = evaluate_fiber(htd, varargin)
%EVALUATE_FIBER Evaluate HTD along a fiber (1 mode vector, others fixed).
%
%   val = ht.evaluate_fiber(htd, id1, id2, ..., idD)
%
%   Input:
%       htd         - HTD structure with fields .U, .B, .tree, .rank
%       id1,...,idD - index arguments
%       Exactly one of the indices must be a vector, the rest must be scalars.
%
%   Output: 
%       val - [N x 1] values of HTD along that fiber.

    if length(varargin) == 1
        % htd is vector
        val = htd.U{varargin{1}};
        return
    end

    tree = htd.tree;
    D    = tree.orders(1);
    if length(varargin) ~= D
        error("Number of input indices must match tensor order %d", D);
    end

    % 1. Identify fiber mode and build idx metadata
    fiber_mode = -1;
    N = -1;
    fixed = true(1, D);
    idx_val = cell(1, D);
    for d = 1:D
        id = varargin{d};
        if ~isscalar(id)
            if fiber_mode ~= -1
                error("Only one input index can be a vector.");
            end
            fiber_mode = d;
            N = length(id);
            fixed(d) = false;
        end
        idx_val{d} = id(:);  % store as column vector
    end
    if fiber_mode == -1
        error("One input must be a vector (fiber mode).");
    end

    % 2. Bottom-up computation
    N_node  = tree.N_node;
    V       = cell(N_node, 1);       % Store values (some scalar, some vector)
    is_fiber= false(N_node, 1);  % Flag whether this node depends on fiber

    % 2.1 Evaluate U at leaf nodes
    for d = 1:D
        node = tree.dim2ind(d);
        U = htd.U{node};  % [N_d x r_d]
        if fixed(d)
            V{node} = U(idx_val{d}, :);  % [1 x r_d]
            is_fiber(node) = false;
        else
            V{node} = U(idx_val{d}, :);  % [N x r_d]
            is_fiber(node) = true;
        end
    end

    % 2.2 Upward pass over internal nodes
    for i = tree.postorder_nonleaf
        l = tree.children(i, 1);
        r = tree.children(i, 2);
        B = htd.B{i};  % [r_l x r_r x r_i]
        r_i = htd.rank(i);
        r_l = htd.rank(l);
        r_r = htd.rank(r);

        Ul = V{l}; Cl = is_fiber(l);
        Ur = V{r}; Cr = is_fiber(r);

        if ~Cl && ~Cr  % Both fixed
            AB          = pagemtimes(Ul,reshape(B,r_l,[]));
            CC          = repmat(Ur, 1, r_i);
            V_          = reshape(AB.*CC,1,r_r,r_i);
            V{i}   = reshape(sum(V_,2),1,r_i);
            is_fiber(i) = false;

        elseif Cl && ~Cr  % Left is fiber
            Q = size(Ul, 1);
            
            if Q*r_r*r_i < 1e5
                % Use vectorized fast path for small size, fallback to loop if too large
                % Threshold ~5e5 elements ≈ 0.8 MB for double precision
                AB   = pagemtimes(Ul,reshape(B,r_l,[]));
                CC   = repmat(Ur, 1, r_i);
                V_   = reshape(AB.*CC,Q,r_r,r_i);
                V{i} = reshape(sum(V_,2),Q,r_i);
            else
                V{i} = zeros(Q, r_i);
                for k = 1:r_i
                    Bk = B(:, :, k);
                    V{i}(:, k) = sum(Ul .* (Ur * Bk.'), 2);
                end
            end
            is_fiber(i) = true;

        elseif ~Cl && Cr  % Right is fiber
            Q = size(Ur, 1); 
            if Q*r_r*r_i < 1e5
                % Use vectorized fast path for small size, fallback to loop if too large
                % Threshold ~1e5 elements ≈ 0.8 MB for double precision
                AB   = pagemtimes(Ul,reshape(B,r_l,[]));
                CC   = repmat(Ur, 1, r_i);
                V_   = reshape(AB.*CC,Q,r_r,r_i);
                V{i} = reshape(sum(V_,2),Q,r_i);
            else
                V{i} = zeros(Q, r_i);
                for k = 1:r_i
                    Bk = B(:, :, k);
                    V{i}(:, k) = sum(Ul .* (Ur * Bk.'), 2);
                end
            end
            is_fiber(i) = true;

        else
            error("Only one branch can depend on fiber — input error?");
        end
    end

    % 3. Final value: V{1} is [N x 1] or [1 x 1]
    val = V{1};
end
