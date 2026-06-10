function htd_new = squeeze(htd)
% SQUEEZE Remove singleton modes via recursive prune-and-merge in one pass.
%
%   htd_new = ht.squeeze(htd)
%
%   This function removes all singleton modes (i.e., dimensions of size 1)
%   from the input HTD by recursively pruning leaf nodes with mode size 1
%   and restructuring the tree. When an interior node is left with only one
%   meaningful child, it is merged into that child. The procedure performs
%   all pruning and merging in a single bottom-up pass.
%
%   Inputs:
%     htd - HTD structure with fields .U, .B, .tree, .rank
%
%   Output:
%     htd_new - squeezed HTD with singleton modes removed
%
%   INTERNAL LOGIC:
%     - Singleton modes are identified by checking mode sizes = 1.
%     - Nodes representing singleton modes are pruned.
%     - Internal nodes with only one meaningful child are merged recursively.
%     - Transfer tensors are updated accordingly.
%
%   NESTED FUNCTION: prune_node(u)
%     Recursively determines which nodes to keep, and performs merging where
%     needed. Uses a reverse node-indexing scheme to reassemble a compact
%     new tree with updated indexing.
%
%   HELPER FUNCTIONS:
%     evaluate_subbasis     - Assembles a new basis via contraction of B, U_l, U_r.
%     mergeNonLeaf_left     - Merges a left child (B_l) into current B.
%     mergeNonLeaf_right    - Merges a right child (B_r) into current B.
%
%   NOTE:
%     This operation is **structure-altering**: it changes the tree topology.
%     It is useful as a post-processing step when some input modes become
%     trivial, e.g., after fixing boundary conditions or solving low-rank PDEs.

% === Step 1: Extract and initialize ===
tree = htd.tree;
N    = tree.N_node;

U_temp        = htd.U;               % working copy of U
B_temp        = htd.B;               % working copy of B
children_temp = tree.children;       % working copy of tree
node_type     = tree.nonleaf;        % node type flag
freedom       = tree.freedom;        % degrees of freedom
isKept        = false(1,N);          % logical: whether node kept
isLeaf        = false(1,N);          % logical: whether node becomes leaf

% === Step 2: Identify modes to keep ===
kept_dims = find(tree.modesizes > 1);   % non-singleton physical dims
L = numel(kept_dims);                   % new number of physical dims
Nnew = 2*L - 1;                         % new number of nodes (full binary tree)

old2new = zeros(N,1);                  % maps old node index to new index
children_new = zeros(Nnew,2);          % new tree structure
freedom_new  = zeros(1,Nnew);          % new freedom vector

next_id = Nnew;                        % assign new node index from bottom up

% === Step 3: Run recursive prune-merge ===
prune_node(1);                         % start from root

% === Step 4: Build new tree structure ===
children = children_new;
old_leaf = tree.dim2ind(kept_dims);
dim2ind  = arrayfun(@(n) old2new(n), old_leaf);
modesizes = tree.modesizes(kept_dims);
tree_new = ht.make_tree(children_new, dim2ind, modesizes);

% === Step 5: Rebuild U and B ===
keepIdx = find(isKept);
U_new = cell(1,Nnew);
B_new = cell(1,Nnew);
rank_new = zeros(1,Nnew);
dim_new2old = find(tree.modesizes > 1);

% Reconstruct U for leaf nodes
for d = 1: tree_new.orders(1)
    node_new = tree_new.dim2ind(d);
    node_old = find(old2new == node_new,1);
    U_new{node_new} = U_temp{node_old};
    rank_new(node_new) = size(U_new{node_new},2);
end

% Reconstruct B for non-leaf nodes
for i = tree_new.postorder_nonleaf
    node_old = find(old2new == i,1);
    B_new{i} = B_temp{node_old};
    rank_new(i) = size(B_new{i},3);
end

% === Step 6: Output new HTD ===
htd_new.tree = tree_new;
htd_new.rank = rank_new;
htd_new.U    = U_new;
htd_new.B    = B_new;


%% ------------------ Nested Recursive Function ------------------ %%
    function [kept, leafFlag, new_idx] = prune_node(u)
        ch = children_temp(u,:);
        % Base Case: original leaf node
        if all(ch==0)
            if freedom(u)>1
                % Keep non-trivial leaf node
                kept      = true;
                leafFlag  = true;
                isKept(u) =true;
                isLeaf(u) =true;
                % U_temp(u) unchanged
                new_idx = next_id;
                old2new(u) = new_idx;
                children_new(new_idx, :) = [0, 0];
                freedom_new(new_idx) = tree.freedom(u);
                next_id = next_id - 1;
            else
                kept = false; leafFlag = false; new_idx = 0;
            end
            return;
        end

        % Recursive pruning: right then left
        [rKept, rLeaf, right_id] = prune_node(ch(2));
        [lKept, lLeaf, left_id] = prune_node(ch(1));

        % Determine which children remain
        rRem = rKept  && freedom(ch(2))>1;
        lRem = lKept  && freedom(ch(1))>1;

        % Default return values
        kept     = false;
        leafFlag = false;
        
        if rRem && lRem
            % Case A: both remain
            kept = true; leafFlag = false;
            isKept(u) = true; isLeaf(u) = false;
            % B_temp(u) unchanged
            new_idx = next_id;
            old2new(u) = new_idx;
            children_new(new_idx, :) = [left_id, right_id];
            freedom_new(new_idx) = tree.freedom(u);
            next_id = next_id - 1;
        elseif xor(rRem,lRem)
            % Case B: only one child kept, try to merge
            kept = true; isKept(u) = true;
            new_idx = max(left_id, right_id); 
            old2new(u) = new_idx;
            if rRem
                keptChild = ch(2);
                childLeaf = rLeaf;
            else
                keptChild = ch(1);
                childLeaf = lLeaf;
            end
            if childLeaf
                % Case B1: merge into new leaf
                leafFlag = true;
                isLeaf(u)=true;
                node_type(u)=false;
                U_temp{u} = evaluate_subbasis(B_temp{u}, ...
                                             U_temp{ch(1)}, ... % left
                                             U_temp{ch(2)});    % right
                isKept(ch(:)) = false;
            else
                % Case B2: merge into new nonleaf
                leafFlag = false;
                isLeaf(u)=false;
                node_type(u)=true;
                % inherit children of keptChild
                children_temp(u,:) = children_temp(keptChild,:);
                % merge B
                if rRem
                    B_temp{u} = mergeNonLeaf_right(B_temp{u}, ...
                                                   U_temp{ch(1)}, ...
                                                   B_temp{ch(2)});
                else
                    B_temp{u} = mergeNonLeaf_left(B_temp{u}, ...
                                                  U_temp{ch(2)}, ...
                                                  B_temp{ch(1)});
                end
                isKept(ch(:)) = false;
            end
        else
            % Case C: none remain
            new_idx = 0;
            if freedom(u)>1
                kept     = true;
                leafFlag = true;
                isKept(u)   =true;
                isLeaf(u)   =true;
                U_temp{u} = evaluate_subbasis(B_temp{u}, ...
                                             U_temp{ch(1)}, ...
                                             U_temp{ch(2)});
            else
                kept = false;
                leafFlag = false;
                isKept(u)   =false;
                isLeaf(u)   =false;
                U_temp{u} = evaluate_subbasis(B_temp{u}, ...
                                             U_temp{ch(1)}, ...
                                             U_temp{ch(2)});
            end
        end
    end

end

%% Helper functions
function U_i = evaluate_subbasis(B, U_l, U_r)
    X = ht.ttm(B, U_l, 1);
    Y = ht.ttm(X, U_r, 2);
    U_i = reshape(Y, [], size(B,3));
end

function B_new = mergeNonLeaf_left(B, U_r, B_l)
    T = ht.ttm(B, U_r, 2);                 % [r_l × 1 × r_i]
    M = reshape(T, size(B,1), []);     % [r_l × r_i]
    B_new = ht.ttm(B_l, M.', 3);          % [r_ll × r_lr × r_i]
end

function B_new = mergeNonLeaf_right(B, U_l, B_r)
    T = ht.ttm(B, U_l, 1);                 % [1 × r_r × r_i]
    M = reshape(T, size(B,2), []);     % [r_r × r_i]
    B_new = ht.ttm(B_r, M.', 3);          % [r_rl × r_rr × r_i]
end
