function htd_new = multiply_scalar(htd,a)
%MULTIPLY_SCALAR Multiply an HTD by a scalar a
%
%   htd_new = ht.multiply_scalar(htd,a)
%
%   Input:
%       htd - HTD structure with fields .U, .B, .tree, .rank
%       a - scalar value to multiply with the HTD
%
%   Output:
%       htd_new - new HTD with root tensor scaled by a

    % Initialize htd_new
    htd_new = htd;

    % Multiply the root transfer tensor of the HTD by scalar a
    htd_new.B{1} = a*htd.B{1};

end