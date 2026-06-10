function total_mass = compute_tmass_2D2V(f,h_mesh)

    B34   = reshape(f.B{3},f.rank(6),[]);
    U34   = sum(f.U{6},1)*B34;
    U4_rp = repmat(sum(f.U{7},1), 1, f.rank(3));
    U34   = U34.*U4_rp;
    U34   = reshape(U34,f.rank(7),f.rank(3));
    U34   = h_mesh(3)*h_mesh(4)*reshape(sum(U34,1),1,f.rank(3));

    B12   = reshape(f.B{2},f.rank(4),[]);
    U12   = sum(f.U{4},1)*B12;
    U2_rp = repmat(sum(f.U{5},1), 1, f.rank(2));
    U12   = U12.*U2_rp;
    U12   = reshape(U12,f.rank(5),f.rank(2));
    U12   = h_mesh(1)*h_mesh(2)*reshape(sum(U12,1),1,f.rank(2));

    total_mass = U12*f.B{1}*U34.';

end