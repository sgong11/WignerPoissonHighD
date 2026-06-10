function [subtree_pool, par2lsub_pool, lsub2par_pool, par2rsub_pool, rsub2par_pool] = get_subtree_pool(tree)
%GET_SUBTREE_POOL Recursively extract all left/right subtrees of a HT tree
%
%   [subtree_pool, par2lsub_pool, lsub2par_pool, par2rsub_pool, rsub2par_pool] = ht.get_subtree_pool(tree)
%
%   Given a full HT tree structure, this function recursively extracts and stores
%   all left and right subtrees rooted at each internal node.
%   The output includes subtree structures and mapping information between
%   each subtree and its parent tree.
%
%   INPUT:
%     tree - full hierarchical Tucker tree (struct with standard fields)
%
%   OUTPUT:
%     subtree_pool  - cell array where subtree_pool{i} is the HT subtree rooted at node i
%     par2lsub_pool - cell array: par2lsub_pool{u} maps original indices in parent to left subtree
%     lsub2par_pool - cell array: lsub2par_pool{i} maps left subtree indices back to parent indices
%     par2rsub_pool - same as above but for right subtree
%     rsub2par_pool - same as above but for right subtree

    N = tree.N_node;            % total number of nodes in the original tree
    children = tree.children;   
    subtree_pool  = cell(1,N);
    lsub2par_pool = cell(1,N);
    rsub2par_pool = cell(1,N);
    par2lsub_pool = cell(1,N);
    par2rsub_pool = cell(1,N);

    % Root node (node 1) is the entire tree
    subtree_pool{1} = tree;

    % Begin recursive extraction of subtrees from root
    get_subtree(1);

    function get_subtree(u)
        ch = children(u,:);
        if all(ch == 0)
            % Base case: leaf node — no subtrees to extract
            return
        else
            % Recursive case: internal node — extract left and right subtrees

            % Extract left subtree rooted at ch(1)
            [subtree_pool{ch(1)} par2lsub_pool{u} lsub2par_pool{ch(1)}] = ht.get_left_subtree(subtree_pool{u});
            % Extract right subtree rooted at ch(2)
            [subtree_pool{ch(2)} par2rsub_pool{u} rsub2par_pool{ch(2)}] = ht.get_right_subtree(subtree_pool{u});

            % Continue recursively on left and right children
            get_subtree(ch(1));
            get_subtree(ch(2));
        end
    end

end