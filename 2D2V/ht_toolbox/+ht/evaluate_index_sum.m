function val_sum = evaluate_index_sum(htd_list, varargin)
%EVALUATE_INDEX_SUM Evaluate and sum multiple HTDs at given multi-indices.
%
%   val_sum = ht.evaluate_index_sum(htd_list, id1, id2, ..., idD)
%
%   Inputs:
%       htd_list - cell array of HTDs {htd1, htd2, ..., htdN}, with identical .tree
%       id1,...,idD - index arguments
%       each of the indix arguments is a [Q x 1] vector containing Q
%       indices of the corresponding mode
%
%   Output:
%       val_sum  - [Q x 1] vector, sum of all evaluations at given indices

    N = numel(htd_list);
    if N == 0
        % The input htd set is empty
        Q = length(varargin{1});
        val_sum  = zeros(Q,1);
        return
    end
    val_sum = zeros(length(varargin{1}), 1); % [Q x 1]

    for k = 1:N
        val = ht.evaluate_index(htd_list{k}, varargin{:}); % [Q x 1]
        val_sum = val_sum + val;
    end
end
