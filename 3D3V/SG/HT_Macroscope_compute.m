function [mass, momentum, energy] = HT_Macroscope_compute(f_ht, h_mesh, x_modes, v_modes, v, E_field,elec_mod)
% HT_MACROSCOPE_COMPUTE Computes macroscopic quantities of HT distribution
%
%   Computes the total mass, momentum vector, and total energy 
%   (kinetic + electric) of the system.

    mass = compute_tmass(f_ht, h_mesh);
    momentum = compute_momentum(f_ht, h_mesh, x_modes, v_modes, v);
    kinetic_energy = real(compute_kinetic_energy(f_ht, h_mesh, x_modes, v_modes, v));
    electric_energy = real(compute_electric_field_energy(h_mesh(x_modes), E_field));
    
    energy = kinetic_energy + electric_energy + elec_mod;
end
