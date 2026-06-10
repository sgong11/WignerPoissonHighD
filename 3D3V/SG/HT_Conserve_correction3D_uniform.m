function [f_ht_new, electric_energy_mod] = HT_Conserve_correction3D_uniform(f_ht, h_mesh, v_modes, x_modes, v, target_mass, target_mom, target_energy, current_mass, current_mom, current_energy, phix_ht_current, base_rel_eps, tol_ratio, min_rank, max_rank, kx)
% HT_CONSERVE_CORRECTION3D_UNIFORM Enforce global conservation on a 3D3V HT tensor
%
%   Applies an additive correction to the tensor f_ht to enforce the 
%   conservation of mass, momentum, and energy in 3D3V.

    % 1. Build the leading-rank approximation used for the correction basis.
    f1_ht = ht.truncate(f_ht,inf,min_rank,min_rank);
    
    
    % 2. Solve for the conservation-correction coefficients.
    [c, eta_tilde_vec, Gamma1_tilde, Gamma2_tilde] = HT_Get_coeff3D(f1_ht, h_mesh, v_modes, x_modes, v, target_mass, target_mom, target_energy, current_mass, current_mom, current_energy, phix_ht_current, base_rel_eps, tol_ratio, min_rank, max_rank, kx);
    
    % 3. Form the correction tensors
    eta = dot(eta_tilde_vec, c);
    electric_energy_mod = eta*Gamma1_tilde + 0.5*eta^2*Gamma2_tilde;
    
    % T1 weights the vz leaf: c1/3 + c4 * vz + c5 * vz^2.
    T1 = f1_ht;
    T1.U{7} = f1_ht.U{7} .* (c(1)/3 + c(4)*v{3}(:) + c(5)*(v{3}(:).^2));
    
    % T2 weights the vx leaf: c1/3 + c2 * vx + c5 * vx^2.
    T2 = f1_ht;
    T2.U{9} = f1_ht.U{9} .* (c(1)/3 + c(2)*v{1}(:) + c(5)*(v{1}(:).^2));
    
    % T3 weights the vy leaf: c1/3 + c3 * vy + c5 * vy^2.
    T3 = f1_ht;
    T3.U{11} = f1_ht.U{11} .* (c(1)/3 + c(3)*v{2}(:) + c(5)*(v{2}(:).^2));
    
    % 4. Add the corrections to f_ht
    f_ht_new = ht.add(f_ht, T1);
    f_ht_new = ht.add(f_ht_new, T2);
    f_ht_new = ht.add(f_ht_new, T3);
end
