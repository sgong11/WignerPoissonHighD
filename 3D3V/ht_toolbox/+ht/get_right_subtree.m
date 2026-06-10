function [tree_r, old2new, new2old] = get_right_subtree(tree)
%GET_RIGHT_SUBTREE Extract the right subtree rooted at tree.children(1,2)
%
%   [tree_l, old2new, new2old] = ht.get_right_subtree(tree)
%
%   This function extracts the right subtree of a given hierarchical Tucker
%   tree structure. The subtree is rooted at the right child of the root
%   node (i.e., tree.children(1,2)). The new subtree is rebuilt with
%   renumbered node indices in a bottom-up fashion.
%
%   INPUT:
%     tree     - original HT tree with fields:
%                .children   : [N×2] connectivity array
%                .dim2ind    : mapping from physical dimensions to leaf nodes
%                .modesizes  : vector of mode sizes
%                .orders     : vector storing the number of modes under each node
%
%   OUTPUT:
%     tree_r   - new tree structure representing the right subtree
%     old2new  - map from original node indices to right-subtree indices (0 if not used)
%     new2old  - inverse map from new subtree indices to original node indices

% === Step 1: Extract basic info from input tree ===
N        = tree.N_node;                  % number of nodes in original tree
children = tree.children;                % original tree structure
root_r   = children(1,2);                % index of right child of root

% === Step 2: Determine right-subtree size and prepare working copies ===
L            = tree.orders(root_r);           % number of physical dimensions in right subtree
Nnew         = 2*L - 1;                       % number of nodes in a full binary tree with L leaves

old2new      = zeros(N,1);                    % map: original index → right-subtree index
new2old      = zeros(Nnew,1);                 % map: right-subtree index → original index
children_new = zeros(Nnew,2);                 % storage for new tree structure

next_id = Nnew;   % node indices assigned in bottom-up order (leaf: Nnew, root: 1)

% === Step 3: Run DFS traversal from root_r to populate new tree ===
recursive_search(root_r);    % populate old2new, new2old, children_new

% === Step 4: Construct new tree using mapped dim2ind and modesizes ===
old_leaf   = tree.dim2ind(end-L+1:end);              % original leaf node indices for dims end-L+1:end
dim2ind    = arrayfun(@(n) old2new(n), old_leaf);   % mapped leaf indices in right subtree
modesizes  = tree.modesizes(end-L+1:end);            % corresponding mode sizes
tree_r     = ht.make_tree(children_new, dim2ind, modesizes);  % assemble output

%% === Nested DFS function for node relabeling ===
    function new_idx = recursive_search(u)
        % Recursively visit subtree rooted at node u and assign new index
        ch = children(u,:);
        if all(ch == 0)
            % Case: leaf node
            new_idx = next_id;
            old2new(u) = new_idx;
            new2old(new_idx) = u;
            children_new(new_idx, :) = [0, 0];
            next_id = next_id - 1;
        else
            % Case: internal node — process children first
            idx_r = recursive_search(ch(2));   % right child
            idx_l = recursive_search(ch(1));   % left child
            new_idx = next_id;
            old2new(u) = new_idx;
            new2old(new_idx) = u;
            children_new(new_idx, :) = [idx_l, idx_r];
            next_id = next_id - 1;
        end
    end

end
