function r = ht_rank(x)
%HT_RANK Compute the rank vector from a given HTD structure.
%
%   r = ht_rank(x)
%
%   Input:
%       x - HTD struct with fields:
%             .U : cell array of leaf basis matrices (size: N_i × r_i)
%             .B : cell array of transfer tensors (size: r_l × r_r × r_i)
%
%   Output:
%       r - [1 x N] row vector, where r(i) is the rank at node i

    N = length(x.U);
    r = zeros(1, N);  % Row vector

    for i = 1:N
        if ~isempty(x.U{i})
            r(i) = size(x.U{i}, 2);  % Leaf node rank
        elseif ~isempty(x.B{i})
            r(i) = size(x.B{i}, 3);  % Internal node rank
        else
            error("Node %d has neither U nor B defined.", i);
        end
    end
end