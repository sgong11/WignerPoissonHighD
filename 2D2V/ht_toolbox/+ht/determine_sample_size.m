function sample_size = determine_sample_size(D)
%DETERMINE_SAMPLE_SIZE Compute the sample size from a given order of tensor.
%
%   sample_size = ht.determine_sample_size(D)
%
%   Input:
%       D - the order of a tensor
%
%   Output:
%       sample_size - scalar, initial sample size used for HTACA
    if D <= 6
        sample_size = 3^D;
    else
        sample_size = 3^6;
    end

end