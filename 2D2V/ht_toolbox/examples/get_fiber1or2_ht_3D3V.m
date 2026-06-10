function f_val = get_fiber1or2_ht_3D3V(f,id_1,id_2,id_3,id_4,id_5,id_6,rk,N)

    B12   = reshape(f.B{4},rk(8),[]);
    U12   = f.U{8}(id_1,:)*B12;
    U2_rp = repmat(f.U{9}(id_2,:), 1, rk(4));
    U12   = U12.*U2_rp;
    U12   = reshape(U12,N,rk(9),rk(4));
    U12   = reshape(sum(U12,2),N,rk(4));

    B34   = reshape(f.B{5},rk(10),[]);
    U34   = f.U{10}(id_3,:)*B34;
    U4_rp = repmat(f.U{11}(id_4,:), 1, rk(5));
    U34   = U34.*U4_rp;
    U34   = reshape(U34,1,rk(11),rk(5));
    U34   = reshape(sum(U34,2),1,rk(5));

    B1234 = reshape(f.B{2},rk(4),[]);
    U1234 = U12*B1234;
    U34_rp= repmat(U34, 1, rk(2));
    U1234 = U1234.*U34_rp;
    U1234 = reshape(U1234,N,rk(5),rk(2));
    U1234 = reshape(sum(U1234,2),N,rk(2));

    B56   = reshape(f.B{3},rk(6),[]);
    U56   = f.U{6}(id_5,:)*B56;
    U6_rp = repmat(f.U{7}(id_6,:), 1, rk(3));
    U56   = U56.*U6_rp;
    U56   = reshape(U56,1,rk(7),rk(3));
    U56   = reshape(sum(U56,2),1,rk(3));
   
    f_val = sum(U1234*f.B{1}.*U56,2);

end