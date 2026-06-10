function f_val = get_val_ht_2D2V(f,id_1,id_2,id_3,id_4,rk)

    N = length(id_1);

    B12   = reshape(f.B{2},rk(4),[]);
    U12   = f.U{4}(id_1,:)*B12;
    U2_rp = repmat(f.U{5}(id_2,:), 1, rk(2));
    U12   = U12.*U2_rp;
    U12   = reshape(U12,N,rk(5),rk(2));
    U12   = reshape(sum(U12,2),N,rk(2));

    B34   = reshape(f.B{3},rk(6),[]);
    U34   = f.U{6}(id_3,:)*B34;
    U4_rp = repmat(f.U{7}(id_4,:), 1, rk(3));
    U34   = U34.*U4_rp;
    U34   = reshape(U34,N,rk(7),rk(3));
    U34   = reshape(sum(U34,2),N,rk(3));
    
    f_val = sum(U12*f.B{1}.*U34,2);

end