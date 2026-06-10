function validate_tree(tree)
%VALIDATE_TREE Checks whether a tree struct is valid and complete.
%
%   ht.validate_tree(tree)
%
%   validate_tree(tree) throws an error if any required field is missing
%   or inconsistent. The expected fields are:
%       - tree.children: [N x 2] array
%       - tree.dim2ind : [1 x D] array of indices into children
%       - tree.modesizes: [1 x D] array of positive integers
%
%   All indices in dim2ind must point to leaf nodes in tree.children.

    % Check required fields
    required_fields = {'children', 'dim2ind', 'modesizes'};
    for f = required_fields
        if ~isfield(tree, f{1})
            error('validate_tree:MissingField', ...
                'Missing required field "%s" in tree structure.', f{1});
        end
    end

    % Check children shape
    children = tree.children;
    if ~ismatrix(children) || size(children, 2) ~= 2
        error('validate_tree:InvalidChildren', ...
            '"children" must be a matrix with two columns.');
    end
    N = size(children, 1);

    % Check dim2ind indices are within range
    if any(tree.dim2ind < 1) || any(tree.dim2ind > N)
        error('validate_tree:InvalidDim2Ind', ...
            '"dim2ind" contains indices outside the valid range of nodes.');
    end

    % Check modesizes length matches dim2ind
    if length(tree.dim2ind) ~= length(tree.modesizes)
        error('validate_tree:MismatchedSizes', ...
            'Length of "dim2ind" must match length of "modesizes".');
    end

    % Check all dim2ind entries point to leaf nodes
    for idx = tree.dim2ind
        if any(children(idx, :) ~= 0)
            error('validate_tree:NonLeafDimIndex', ...
                'dim2ind maps to non-leaf node at index %d.', idx);
        end
    end

    % Check modesizes entries are positive integers
    if any(tree.modesizes <= 0) || any(mod(tree.modesizes, 1) ~= 0)
        error('validate_tree:InvalidModesizes', ...
            '"modesizes" must contain positive integers.');
    end
end
