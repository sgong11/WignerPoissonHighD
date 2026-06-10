function phi = eval_dir_clamped_3D(signH, k1_vec, k2_vec, k3_vec, dx, dy, dz, ...
                                   idx_x, idx_y, idx_z, ...
                                   phi_ht, phix_ht, phi_cross, x_extend, n_ghost, N_mesh)
    % 1. Calculate shifted interpolation locations.
    p_x = idx_x - signH * k1_vec / dx;
    p_y = idx_y - signH * k2_vec / dy;
    p_z = idx_z - signH * k3_vec / dz;

    % 2. Base cell node.
    ix_float = floor(p_x);
    iy_float = floor(p_y);
    iz_float = floor(p_z);

    % 3. Fractional distances [0, 1)
    t = p_x - ix_float;
    u = p_y - iy_float;
    v = p_z - iz_float;

    % 4. Periodic wrap for the base node
    ix = mod(ix_float - 1, N_mesh(1)) + 1;
    iy = mod(iy_float - 1, N_mesh(3)) + 1;
    iz = mod(iz_float - 1, N_mesh(5)) + 1;

    % 5. 1D Hermite basis functions on [0, 1].
    h00_t = 2*t.^3 - 3*t.^2 + 1;  h10_t = t.^3 - 2*t.^2 + t;
    h01_t = -2*t.^3 + 3*t.^2;     h11_t = t.^3 - t.^2;

    h00_u = 2*u.^3 - 3*u.^2 + 1;  h10_u = u.^3 - 2*u.^2 + u;
    h01_u = -2*u.^3 + 3*u.^2;     h11_u = u.^3 - u.^2;

    h00_v = 2*v.^3 - 3*v.^2 + 1;  h10_v = v.^3 - 2*v.^2 + v;
    h01_v = -2*v.^3 + 3*v.^2;     h11_v = v.^3 - v.^2;

    H_t = {h00_t, h01_t}; G_t = {h10_t, h11_t};
    H_u = {h00_u, h01_u}; G_u = {h10_u, h11_u};
    H_v = {h00_v, h01_v}; G_v = {h10_v, h11_v};

    phi = zeros(size(p_x));

    % 6. Evaluate over the 8 corners directly using spectral derivatives
    for cx = 0:1
        for cy = 0:1
            for cz = 0:1
                % Periodic corner indices.
                curr_ix = mod(ix + cx - 1, N_mesh(1)) + 1;
                curr_iy = mod(iy + cy - 1, N_mesh(3)) + 1;
                curr_iz = mod(iz + cz - 1, N_mesh(5)) + 1;

                eval_x = x_extend{1}(n_ghost + curr_ix, 2);
                eval_y = x_extend{2}(n_ghost + curr_iy, 2);
                eval_z = x_extend{3}(n_ghost + curr_iz, 2);

                % Batch tensor evaluation for N coordinates at this specific cell corner
                val_phi = real(ht.evaluate_index(phi_ht, eval_x, eval_y, eval_z));
                val_x   = real(ht.evaluate_index(phix_ht{1}, eval_x, eval_y, eval_z));
                val_y   = real(ht.evaluate_index(phix_ht{2}, eval_x, eval_y, eval_z));
                val_z   = real(ht.evaluate_index(phix_ht{3}, eval_x, eval_y, eval_z));
                val_xy  = real(ht.evaluate_index(phi_cross.xy, eval_x, eval_y, eval_z));
                val_xz  = real(ht.evaluate_index(phi_cross.xz, eval_x, eval_y, eval_z));
                val_yz  = real(ht.evaluate_index(phi_cross.yz, eval_x, eval_y, eval_z));
                val_xyz = real(ht.evaluate_index(phi_cross.xyz, eval_x, eval_y, eval_z));

                % Map derivative terms back from cell-normalized coordinates.
                wx = H_t{cx+1}; wgx = G_t{cx+1} * dx;
                wy = H_u{cy+1}; wgy = G_u{cy+1} * dy;
                wz = H_v{cz+1}; wgz = G_v{cz+1} * dz;

                % Add this corner contribution.
                term = val_phi .* (wx .* wy .* wz) + ...
                       val_x   .* (wgx .* wy .* wz) + ...
                       val_y   .* (wx .* wgy .* wz) + ...
                       val_z   .* (wx .* wy .* wgz) + ...
                       val_xy  .* (wgx .* wgy .* wz) + ...
                       val_xz  .* (wgx .* wy .* wgz) + ...
                       val_yz  .* (wx .* wgy .* wgz) + ...
                       val_xyz .* (wgx .* wgy .* wgz);

                phi = phi + term;
            end
        end
    end
end
