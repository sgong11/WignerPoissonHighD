function phihalf = get_phi_clamp(phi_ht,phix_ht,phixy_ht,l,index,x_extend, n_ghost, ...
                                 k_mod_1,k_mod_2,H,h_mesh, ...
                                 N_mesh,dom_len)

    dx = h_mesh(1);
    dy = h_mesh(3);

    idx = n_ghost+index(l,:);

    
    % Scaling for shifts
    theta_x = (H / 2) * k_mod_1(:)'; 
    theta_y = (H / 2) * k_mod_2(:)';

    xp = x_extend{1}(idx(1)) + theta_x;
    xm = x_extend{1}(idx(1)) - theta_x;
    yp = x_extend{2}(idx(3)) + theta_y;
    ym = x_extend{2}(idx(3)) - theta_y;

    xp_mod = mod(xp - x_extend{1}(n_ghost+1), dom_len(1) ) + x_extend{1}(n_ghost+1);
    xm_mod = mod(xm - x_extend{1}(n_ghost+1), dom_len(1) ) + x_extend{1}(n_ghost+1);
    yp_mod = mod(yp - x_extend{2}(n_ghost+1), dom_len(3) ) + x_extend{2}(n_ghost+1);
    ym_mod = mod(ym - x_extend{2}(n_ghost+1), dom_len(3) ) + x_extend{2}(n_ghost+1);

    PHI_P = interpolate2D_clamped(xp_mod, yp_mod, x_extend, n_ghost, phi_ht, phix_ht,phixy_ht, dx, dy, N_mesh);
    PHI_M = interpolate2D_clamped(xm_mod, ym_mod, x_extend, n_ghost, phi_ht, phix_ht, phixy_ht,dx, dy, N_mesh);
    
    phihalf = PHI_P - PHI_M;
end