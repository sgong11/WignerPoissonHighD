function f_new = SL_HT2D2V_handle_vec(i,j,p,q,f_preT, ...
        x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt,epsilon)

    % normalize inputs
    vars = {i,j,p,q};
    vars = cellfun(@(x) x(:), vars, 'UniformOutput', false);
    [i,j,p,q] = deal(vars{:});
    N = max(cellfun(@length, vars));

    index = [i.*ones(N,1), j.*ones(N,1), p.*ones(N,1), q.*ones(N,1)];

    % coordinates
    coor = zeros(N,4);
    for d = 1:2
        id = x_modes(d);
        coor(:,id) = x_extend{d}(n_ghost + index(:,id), 1);
    end
    for d = 1:2
        id = v_modes(d);
        coor(:,id) = v{d}(index(:,id));
    end

    % precompute shifts (NO while-loops)
    xshift = coor(:,v_modes(1)) * dt / h_mesh(x_modes(1));
    yshift = coor(:,v_modes(2)) * dt / h_mesh(x_modes(2));

    sx = round(xshift);
    sy = round(yshift);

    xshift = xshift - sx;
    yshift = yshift - sy;

    xshift_pos = abs(xshift);
    yshift_pos = abs(yshift);

    f_new = zeros(N,1);

    for l = 1:N
        % shifted indices
        inx_mod = index(l,x_modes(1)) - sx(l);
        iny_mod = index(l,x_modes(2)) - sy(l);

        % x-stencil
        if xshift(l) >= 0
            stencil_x = inx_mod-3:inx_mod+2;
            caseinx = 1;
        else
            stencil_x = inx_mod-2:inx_mod+3;
            caseinx = 2;
        end

        % y-stencil
        if yshift(l) >= 0
            stencil_y = iny_mod-3:iny_mod+2;
            caseiny = 1;
        else
            stencil_y = iny_mod-2:iny_mod+3;
            caseiny = 2;
        end

        % build stencil grid
        [ix, iy] = ndgrid(stencil_x, stencil_y);
        Ns = numel(ix);

        % vectorized HT evaluation
        f_vals = get_val_ht_2D2V_SG( ...
            f_preT, ...
            x_extend{1}(n_ghost + ix(:),2), ...
            x_extend{2}(n_ghost + iy(:),2), ...
            repmat(index(l,v_modes(1)), Ns, 1), ...
            repmat(index(l,v_modes(2)), Ns, 1) );

        f_temp = reshape(real(f_vals), 6, 6);

        % WENO reconstruction
        f_temp_x = SL1DlocalWENO_vec( ...
            f_temp, xshift(l), xshift_pos(l), epsilon, caseinx);
        
        f_new(l) = SL1DlocalWENO_vec( ...
            f_temp_x.', yshift(l), yshift_pos(l), epsilon, caseiny);

    end
end
