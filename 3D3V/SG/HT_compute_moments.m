function M_pqr = HT_compute_moments(f_ht, h_mesh, v_modes, v, p, q, r)
% HT_COMPUTE_MOMENTS Calculates the generalized moments of an HT tensor (3D3V)
% 
%   M_pqr = HT_compute_moments(f_ht, h_mesh, v_modes, v, p, q, r)
%   
%   Computes M_{p,q,r} = \int \int f_ht * v_x^p * v_y^q * v_z^r dx dv
%
%   Inputs:
%       f_ht    - HT tensor
%       h_mesh  - grid spacing
%       v_modes - velocity mode indices (expected length 3)
%       v       - velocity grid (cell array, v{1} = vx, v{2} = vy, v{3} = vz)
%       p       - power of vx
%       q       - power of vy
%       r       - power of vz

    D = f_ht.tree.orders(1);
    
    F_list = cell(1, D);
    index_ = cell(1, D);
    for d = 1:D
        F_list{d}  = @(x) h_mesh(d) * sum(x, 1);
        index_{d}  = 1;
    end
    
    if p > 0
        F_list{v_modes(1)} = @(x) h_mesh(v_modes(1)) * sum(x .* (v{1}(:).^p), 1);
    end
    
    if q > 0
        F_list{v_modes(2)} = @(x) h_mesh(v_modes(2)) * sum(x .* (v{2}(:).^q), 1);
    end
    
    if r > 0
        F_list{v_modes(3)} = @(x) h_mesh(v_modes(3)) * sum(x .* (v{3}(:).^r), 1);
    end
    
    M_pqr = ht.evaluate_index(ht.mode_apply_function(f_ht, F_list, 1:D), index_{:});
end