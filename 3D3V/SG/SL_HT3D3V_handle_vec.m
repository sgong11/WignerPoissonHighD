function f_new = SL_HT3D3V_handle_vec(i,j,p,q,k,r,f_preT, ...
        x_extend,v,h_mesh,x_modes,v_modes,n_ghost,dt,epsilon)

    % Normalize all input indices to column vectors.
    vars = {i,j,p,q,k,r};
    vars = cellfun(@(x) x(:), vars, 'UniformOutput', false);
    [i,j,p,q,k,r] = deal(vars{:});
    N = max(cellfun(@length, vars));

    % index structure: [x, vx, y, vy, z, vz]
    index = [i.*ones(N,1), j.*ones(N,1), ...
             p.*ones(N,1), q.*ones(N,1), ...
             k.*ones(N,1), r.*ones(N,1)];

    % coordinates
    coor = zeros(N,6);
    % Spatial coordinates (x, y, z)
    for d = 1:3
        id = x_modes(d);
        coor(:,id) = x_extend{d}(n_ghost + index(:,id), 1);
    end
    % Velocity coordinates (vx, vy, vz)
    for d = 1:3
        id = v_modes(d);
        coor(:,id) = v{d}(index(:,id));
    end

    % Precompute spatial shifts in x, y, and z.
    xshift = coor(:,v_modes(1)) * dt / h_mesh(x_modes(1));
    yshift = coor(:,v_modes(2)) * dt / h_mesh(x_modes(2));
    zshift = coor(:,v_modes(3)) * dt / h_mesh(x_modes(3));

    sx = round(xshift);
    sy = round(yshift);
    sz = round(zshift);

    xshift = xshift - sx;
    yshift = yshift - sy;
    zshift = zshift - sz;

    xshift_pos = abs(xshift);
    yshift_pos = abs(yshift);
    zshift_pos = abs(zshift);

    f_new = zeros(N,1);

    for l = 1:N
        % Shifted indices at the characteristic foot.
        inx_mod = index(l,x_modes(1)) - sx(l);
        iny_mod = index(l,x_modes(2)) - sy(l);
        inz_mod = index(l,x_modes(3)) - sz(l);

        % x-stencil
        if xshift(l) >= 0
            stencil_x = inx_mod-3:inx_mod+2; caseinx = 1;
        else
            stencil_x = inx_mod-2:inx_mod+3; caseinx = 2;
        end

        % y-stencil
        if yshift(l) >= 0
            stencil_y = iny_mod-3:iny_mod+2; caseiny = 1;
        else
            stencil_y = iny_mod-2:iny_mod+3; caseiny = 2;
        end

        % z-stencil
        if zshift(l) >= 0
            stencil_z = inz_mod-3:inz_mod+2; caseinz = 1;
        else
            stencil_z = inz_mod-2:inz_mod+3; caseinz = 2;
        end

        % Build the 3D interpolation stencil.
        [ix, iy, iz] = ndgrid(stencil_x, stencil_y, stencil_z);
        Ns = numel(ix);

        % Vectorized HT evaluation at all stencil points.
        f_vals = ht.evaluate_index( ...
            f_preT, ...
            x_extend{1}(n_ghost + ix(:),2), ...
            repmat(index(l,v_modes(1)), Ns, 1), ...
            x_extend{2}(n_ghost + iy(:),2), ...
            repmat(index(l,v_modes(2)), Ns, 1), ...
            x_extend{3}(n_ghost + iz(:),2), ...
            repmat(index(l,v_modes(3)), Ns, 1) );

        % Reshape to 6x6x6 cube
        f_temp = reshape(real(f_vals), 6, 6, 6);

        % WENO Reconstruction Sequence:
        % 1. Contract X dimension (6x6x6 -> 1x6x6)
        % Flatten y-z dims to vectorize the x-interp
        f_temp_x = SL1DlocalWENO_vec( ...
            reshape(f_temp, 6, 36), xshift(l), xshift_pos(l), epsilon, caseinx);
        
        % 2. Contract Y dimension (1x6x6 -> 6x6)
        % Reshape result to (6,6) where rows are Y
        f_temp_xy = reshape(f_temp_x, 6, 6); 
        f_temp_y = SL1DlocalWENO_vec( ...
            f_temp_xy, yshift(l), yshift_pos(l), epsilon, caseiny);
        
        % 3. Contract Z dimension (1x6 -> 6x1)
        % Result is now a vector along Z
        f_temp_y = f_temp_y(:); 
        f_new(l) = SL1DlocalWENO_vec( ...
            f_temp_y, zshift(l), zshift_pos(l), epsilon, caseinz);

    end
end
