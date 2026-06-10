function r = choose_rank_fro(s, tol, min_rank, max_rank)
%CHOOSE_RANK_FRO Truncates singular values by Frobenius norm threshold.
%
%   r = ht.choose_rank_fro(S, tol, max_rank)
%
%   Inputs:
%       s        - Vector of singular values (e.g., from svd)
%       tol      - Frobenius norm tolerance (absolute, not relative)
%       min_rank - Lower bound on the allowed rank
%       max_rank - Upper bound on the allowed rank
%
%   Output:
%       r - Truncation rank such that the discarded tail of singular values
%           has Frobenius norm ≤ tol, and r ≤ max_rank. If this is not 
%           possible, it uses max_rank and prints a warning.

    % Enforce rank constraint (initial)
    if length(s) < min_rank
        r = length(s);
        return
    end

    % Compute cumulative Frobenius norm from tail of singular values
    % Flip and accumulate squared singular values from smallest to largest
    s_sum = sqrt(cumsum(s(end:-1:1).^2));
    s_sum = s_sum(end:-1:1);  % Flip back to match original order

    % Find the smallest r such that tail norm < tol
    r = find(s_sum < tol, 1, 'first');

    if isempty(r)
        % No truncation satisfies the tolerance: use full rank
        r = length(s);
    elseif r == 1
        % Keep at least one singular value
        r = 1;
    else
        % Retain singular values up to the one before the violation
        r = r - 1;
    end

    % Enforce rank constraint
    if r > max_rank
        r = max_rank;
        % disp('Tolerance is not satisfied'); % SG deleted
    elseif r < min_rank && min_rank <= length(s)
        r = min_rank;
    end
end
