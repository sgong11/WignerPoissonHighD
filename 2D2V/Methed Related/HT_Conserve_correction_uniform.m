function [f_ht_new, electric_energy_mod] = HT_Conserve_correction_uniform(f_ht, h_mesh, v_modes, x_modes, v, target_mass, target_mom, target_energy, current_mass, current_mom, current_energy,min_rank)
% HT_CONSERVE_CORRECTION_UNIFORM Enforces global conservation using h(x)=1
%
%   Applies a spatially uniform additive correction to the tensor f_ht.
%   This completely avoids solving the Poisson equation for the correction step.
    f1_ht = ht.truncate(f_ht,inf,min_rank,min_rank);
    % 1. Get the purely linear coefficients
    c = HT_Get_coeff_uniform(f1_ht, h_mesh, v_modes, v, target_mass, target_mom(1), target_mom(2), target_energy, current_mass, current_mom(1), current_mom(2), current_energy);
    
    % 2. Electric energy modification is exactly zero
    electric_energy_mod = 0;
    
    % 3. Create the base uniform tensor (rank-1, all ones)
    % T_ones = ht.ones(f_ht.tree);
    
    % 4. Form the correction tensors
    % Note: Using U{5} and U{7} to match the velocity leaves in your tree setup
    T1 = f1_ht;
    T1.U{5} = T1.U{5} .* (c(1)/2 + c(2)*v{1}(:) + c(4)*(v{1}(:).^2));
    
    T2 = f1_ht;
    T2.U{7} = T2.U{7} .* (c(1)/2 + c(3)*v{2}(:) + c(4)*(v{2}(:).^2));
    
    % 5. Add the uniform corrections to f_ht
    f_ht_new = ht.add(f_ht, T1);
    f_ht_new = ht.add(f_ht_new, T2);
end