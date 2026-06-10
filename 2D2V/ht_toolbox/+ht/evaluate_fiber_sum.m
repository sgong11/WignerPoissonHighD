function val_sum = evaluate_fiber_sum(htd_list, varargin)
%EVALUATE_FIBER_SUM Evaluate and sum multiple HTDs along a given fiber.
%
%   val_sum = ht.evaluate_fiber_sum(htd_list, id1, id2, ..., idD)
%
%   Inputs:
%       htd_list - cell array of HTDs {htd1, htd2, ..., htdN}, with identical .tree
%       id1,...,idD - index arguments (same format as evaluate_fiber)
%
%   Output:
%       val_sum - sum of all fiber evaluations

    N = numel(htd_list);
    if N == 0
        % The input htd set is empty
        size_val = 0;
        for k = 1: length(varargin)
            size_val = max(size_val,length(varargin{k}));
        end
        val_sum = zeros(size_val,1);
        return
    end
    % Determine fiber output size
    val0 = ht.evaluate_fiber(htd_list{1}, varargin{:});
    val_sum = zeros(length(val0),1);

    val_sum = val_sum + val0;  % add first one
    for k = 2:N
        val = ht.evaluate_fiber(htd_list{k}, varargin{:});
        val_sum = val_sum + val;
    end
end
