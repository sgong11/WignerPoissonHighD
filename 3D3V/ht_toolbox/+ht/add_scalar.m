function htd_new = add_scalar(htd,a)
%ADD_SCALAR Add a scalar a to all entries of an HTD
%
%   htd_new = ht.add_scalar(htd, a)
%
%   This function adds a scalar value a to every element of an HTD.
%
%   INPUT:
%       htd - HTD structure with fields .U, .B, .tree, .rank
%       a   - scalar value to add (broadcasted across entire tensor)
%
%   OUTPUT:
%       htd_new - HTD structure representing (htd + a)
%
%   NOTE:
%       This is implemented as: htd + a * ht.ones(htd.tree), where
%       ht.ones(htd.tree) constructs an all-one tensor in HTD format.

    % Add scalar 'a' by adding scaled all-ones HTD tensor
    htd_new = ht.add(htd,ht.multiply_scalar(ht.ones(htd.tree),a));

end