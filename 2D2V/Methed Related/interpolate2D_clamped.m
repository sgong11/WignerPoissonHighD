% function val = interpolate2D_clamped(xq, yq, x, y, phi, phi_sx, phi_sy, dx, dy, Nx, Ny)
function val = interpolate2D_clamped(xq, yq, x_extend, n_ghost, phi_ht, phix_ht,phixy_ht, dx, dy, N_mesh)
  
    % Find grid indices
    ix1 = floor((xq - x_extend{1}(n_ghost+1)) / dx) + 1;
    iy1 = floor((yq - x_extend{2}(n_ghost+1)) / dy) + 1;
    ix2 = ix1 + 1;
    iy2 = iy1 + 1;

    % Periodic index mapping
    ixW = [mod(ix1-1, N_mesh(1)) + 1; mod(ix2-1, N_mesh(1)) + 1];
    iyW = [mod(iy1-1, N_mesh(3)) + 1; mod(iy2-1, N_mesh(3)) + 1];

    % Local coordinates [0, 1]
    tx = (xq - x_extend{1}(n_ghost+ix1)) / dx;
    ty = (yq - x_extend{2}(n_ghost+iy1)) / dy;

     % 1D Hermite Basis Functions
    h0 = @(t) 2*t.^3 - 3*t.^2 + 1;
    h1 = @(t) -2*t.^3 + 3*t.^2;
    g0 = @(t) t.^3 - 2*t.^2 + t;
    g1 = @(t) t.^3 - t.^2;

    % Evaluate basis functions for both dimensions
    HX = {h0(tx), h1(tx)}; GX = {g0(tx), g1(tx)};
    HY = {h0(ty), h1(ty)}; GY = {g0(ty), g1(ty)};

    % Initialize result
    val = zeros(size(xq));

    % Symmetric Double Summation (4 corners x 4 types of data = 16 terms)
    for i = 1:2 % x-direction corners
        for j = 1:2 % y-direction corners
            % Corner indices
            row = ixW(i, :); 
            col = iyW(j, :);
            
            % Add Function Value Contribution
            val = val + HX{i} .* HY{j} .* ht.evaluate_index(phi_ht,row, col);
            
            % Add X-Slope Contribution (scaled by dx)
            val = val + GX{i} .* HY{j} .* (dx * ht.evaluate_index(phix_ht{1},row, col));
            
            % Add Y-Slope Contribution (scaled by dy)
            val = val + HX{i} .* GY{j} .* (dy * ht.evaluate_index(phix_ht{2},row, col));
            
            % Add Twist Contribution (scaled by dx*dy)
            val = val + GX{i} .* GY{j} .* (dx * dy * ht.evaluate_index(phixy_ht,row, col));
        end
    end
    
    % phi1 = ht.evaluate_index(phi_ht,ix1w,iy1w); % phi(ix1w, iy1w)
    % phi2 = ht.evaluate_index(phi_ht,ix2w,iy1w); % phi(ix2w, iy1w)
    % phi3 = ht.evaluate_index(phi_ht,ix1w,iy2w); % phi(ix1w, iy2w)
    % phi4 = ht.evaluate_index(phi_ht,ix2w,iy2w); % phi(ix2w, iy2w)
    % phix1 = ht.evaluate_index(phix_ht,ix1w,iy1w); % phix(ix1w, iy1w)
    % phix2 = ht.evaluate_index(phix_ht,ix2w,iy1w); % phix(ix2w, iy1w)
    % phix3 = ht.evaluate_index(phix_ht,ix1w,iy2w); % phix(ix1w, iy2w)
    % phix4 = ht.evaluate_index(phix_ht,ix2w,iy2w); % phix(ix2w, iy2w)
    % 
    % 
    % % Bicubic Interpolation (tensor product of 1D Hermite)
    % % val = sum_{i,j=0}^1 [ phi * H_i(u)H_j(v) + phi_sx * dx * G_i(u)H_j(v) + phi_sy * dy * H_i(u)G_j(v) ]
    % % Indexing: (ix_w, iy_w)
    % val = hu0.*hv0.*phi1 + hu1.*hv0.*phi2 + ...
    %       hu0.*hv1.*phi3 + hu1.*hv1.*phi4 + ...
    %       (dx*phi_sx(ix1w,iy1w)).*gu0.*hv0 + (dx*phi_sx(ix2w,iy1w)).*gu1.*hv0 + ...
    %       (dy*phi_sy(ix1w,iy1w)).*hu0.*gv0 + (dy*phi_sy(ix1w,iy2w)).*hu0.*gv1; 
          % Note: For full accuracy, add the remaining 8 terms for sx and sy at all corners
end