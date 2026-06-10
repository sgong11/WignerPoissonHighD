function [f_ht_new,electric_energy_mod] = HT_Conserve_correction_leading(f_ht, h_mesh, v_modes,x_modes, v, target_mass, target_mom, target_energy, current_mass, current_mom, current_energy, phix_ht_current,base_rel_eps,tol_ratio, min_rank, max_rank,kx)
% HT_CONSERVE_CORRECTION Enforces global conservation on the HT tensor
%
%   Applies an additive correction to the tensor f_ht to enforce the 
%   conservation of mass, momentum, and energy.

    % 1. Extract the "first mode" f1_ht from the root tensor
    f1_ht = ht.truncate(f_ht,inf,min_rank,min_rank);
    % f1_ht.U{5}(:,2:end) = 0;
    % f1_ht.U{7}(:,2:end) = 0;
    
    % 2. Get the coefficients
    [c,eta_tilde_vec,Gamma1_tilde,Gamma2_tilde] = HT_Get_coeff2(f1_ht, h_mesh, v_modes,x_modes, v, target_mass, target_mom(1), target_mom(2), target_energy, current_mass, current_mom(1), current_mom(2), current_energy, phix_ht_current, base_rel_eps, tol_ratio, min_rank, max_rank, kx);
    % 3. Form the correction tensors
    eta = dot(eta_tilde_vec, c);
    electric_energy_mod = eta*Gamma1_tilde + 0.5*eta^2*Gamma2_tilde;
    % T1 = c1 * f1 * vx
    T1 = f1_ht;
    T1.U{5} = f1_ht.U{5} .*(c(1)/2 + c(2)*v{1}(:) + c(4)*(v{1}(:).^2));
    
    
    T2 = f1_ht;
    T2.U{7} = f1_ht.U{7} .*(c(1)/2 + c(3)*v{2}(:) + c(4)*(v{2}(:).^2));
    
    
    % 4. Add the corrections to f_ht
    f_ht_new = ht.add(f_ht, T1);
    f_ht_new = ht.add(f_ht_new, T2);
end
