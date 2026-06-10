function A = evaluate_slice(htd,varargin)
%EVALUATE_SLICE Evaluate a 2D slice of an HTD tensor.
%
%   A = ht.evaluate_slice(htd, id1, id2, ..., idD)
%
%   Given a hierarchical Tucker tensor (HTD format), this function evaluates
%   a 2D slice by varying two specified modes (as vectors) and fixing all other
%   modes (as scalars).
%
%   Input:
%       htd         - HTD structure with fields .U, .B, .tree, .rank
%       id1,...,idD - a total of D inputs (matching the tensor order). Exactly two
%                    of them must be non-scalar vectors (varying slice directions),
%                    and the rest must be scalars (fixed indices).
%
%   Output: 
%       val - a [m × n] matrix containing the evaluated 2D slice of the tensor
%           corresponding to the two varying modes.
%
%   Notes:
%       - If the tensor is 2D (i.e., a matrix), the full matrix is returned.
%       - Errors are raised if the number of inputs does not match the tensor
%         order or if the number of vector inputs is not exactly two.

    % Special case when htd is a 2nd-order tensor (matrix)
    if length(varargin) == 2
        % htd is a matrix
        dim2ind = htd.tree.dim2ind;
        A       = htd.U{dim2ind(1)}*htd.B{1}*(htd.U{dim2ind(2)}).';
        return
    end

    tree = htd.tree;
    D    = tree.orders(1);
    if length(varargin) ~= D
        error("Number of input indices must match tensor order %d", D);
    end

    % Identify slice modes
    idx_vecs = find_two_vectors(varargin);

    % Determine the size of the slice
    m = length(varargin{idx_vecs(1)});
    n = length(varargin{idx_vecs(2)});

    % Assemble the slice
    A     = zeros(m,n);
    index = varargin; 
    for k = 1: n
        index{idx_vecs(2)} = varargin{idx_vecs(2)}(k);

        A(:,k) = ht.evaluate_fiber(htd,index{:});
    end

end

function [idx_vecs] = find_two_vectors(cell_array)
%FIND_TWO_VECTORS Find positions of two vector entries in a cell array.
%
%   idx_vecs = find_two_vectors(cell_array)
%
%   Input:
%       cell_array - 1 x D or D x 1 cell array containing scalars and exactly two vector arrays
%
%   Output:
%       idx_vecs   - 1 x 2 vector of indices where the vector arrays are located
%
%   Errors:
%       - If there are more or fewer than two vector arrays, an error is thrown.

    % Check for vector arrays
    is_vector = cellfun(@(x) isvector(x) && ~isscalar(x), cell_array);

    % Find indices of vector arrays
    idx_vecs = find(is_vector);

    % Error checks
    if numel(idx_vecs) < 2
        error('Expected exactly two vector arrays, but found less than two.');
    elseif numel(idx_vecs) > 2
        error('Expected exactly two vector arrays, but found more than two.');
    end
end