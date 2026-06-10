function [phi, phix_ht, phi_cross] = get_phi_FFT(rho, base_rel_tol, tol_ratio, min_rank, max_rank, kx)
    % Enforce the solvability condition.
    rho = ht.add_scalar(rho, -1);
    
    % Apply FFT to the bases of \rho
    tree = rho.tree;
    D    = tree.orders(1);
    for d = 1:D
        node = tree.dim2ind(d);
        rho.U{node} = fft(rho.U{node}, tree.modesizes(d), 1);
    end
    
    % Construct the function handle for solving Δφ = ρ in the frequency domain
    if D == 2
        data = @(i,j) phi_hat_fnc(rho, kx, D, i, j); 
    else
        data = @(i,j,k) phi_hat_fnc(rho, kx, D, i, j, k); 
    end

    % Construct φ̂ in HT format using HTACA
    phi_hat = ht.HTACA(data, tree, base_rel_tol, tol_ratio, min_rank, max_rank);
    
    % Apply inverse FFT to recover φ
    phi = phi_hat;
    for d = 1:D
        node = tree.dim2ind(d);
        phi.U{node} = ifft(phi_hat.U{node}, tree.modesizes(d), 1);
    end

    % Compute first derivatives of the potential.
    phix_ht = cell(1,D);
    phix_ht = cellfun(@(x) phi, phix_ht, 'UniformOutput', false);
    for d = 1:D
        node = tree.dim2ind(d);
        phix_ht{d}.U{node} = ifft(phi_hat.U{node} .* (1i*kx{d}), tree.modesizes(d), 1);
    end
    
    % Compute exact spectral cross-derivatives for 3D Hermite interpolation.
    if D >= 3
        node_x = tree.dim2ind(1);
        node_y = tree.dim2ind(2);
        node_z = tree.dim2ind(3);

        % XY
        phi_cross.xy = phi_hat;
        phi_cross.xy.U{node_x} = ifft(phi_hat.U{node_x} .* (1i*kx{1}), tree.modesizes(1), 1);
        phi_cross.xy.U{node_y} = ifft(phi_hat.U{node_y} .* (1i*kx{2}), tree.modesizes(2), 1);
        phi_cross.xy.U{node_z} = ifft(phi_hat.U{node_z}, tree.modesizes(3), 1);

        % XZ
        phi_cross.xz = phi_hat;
        phi_cross.xz.U{node_x} = ifft(phi_hat.U{node_x} .* (1i*kx{1}), tree.modesizes(1), 1);
        phi_cross.xz.U{node_y} = ifft(phi_hat.U{node_y}, tree.modesizes(2), 1);
        phi_cross.xz.U{node_z} = ifft(phi_hat.U{node_z} .* (1i*kx{3}), tree.modesizes(3), 1);

        % YZ
        phi_cross.yz = phi_hat;
        phi_cross.yz.U{node_x} = ifft(phi_hat.U{node_x}, tree.modesizes(1), 1);
        phi_cross.yz.U{node_y} = ifft(phi_hat.U{node_y} .* (1i*kx{2}), tree.modesizes(2), 1);
        phi_cross.yz.U{node_z} = ifft(phi_hat.U{node_z} .* (1i*kx{3}), tree.modesizes(3), 1);

        % XYZ
        phi_cross.xyz = phi_hat;
        phi_cross.xyz.U{node_x} = ifft(phi_hat.U{node_x} .* (1i*kx{1}), tree.modesizes(1), 1);
        phi_cross.xyz.U{node_y} = ifft(phi_hat.U{node_y} .* (1i*kx{2}), tree.modesizes(2), 1);
        phi_cross.xyz.U{node_z} = ifft(phi_hat.U{node_z} .* (1i*kx{3}), tree.modesizes(3), 1);
    else
        phi_cross = [];
    end
end

function value = phi_hat_fnc(rho, k_fre, D, varargin)
    index_len = [];
    K = 0;
    for d = 1:D
        index_len = [index_len length(varargin{d})];
        K = K + (k_fre{d}(varargin{d})).^2;
    end
    N = max(index_len);

    one_over_K = zeros(N,1);
    one_over_K(K>0) = 1./K(K>0);

    if sum(index_len > 1) == 1
        value = one_over_K .* ht.evaluate_fiber(rho, varargin{:});
        return
    else
        varargin = cellfun(@(x) x.*ones(N,1), varargin, 'UniformOutput', false);
        value = one_over_K .* ht.evaluate_index(rho, varargin{:});
    end
end
