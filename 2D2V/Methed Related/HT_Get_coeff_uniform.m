function [c] = HT_Get_coeff_uniform(f_ht, h_mesh, v_modes, v, target_mass, target_mom_x, target_mom_y, target_energy, current_mass, current_mom_x, current_mom_y, current_energy)
% HT_GET_COEFF_UNIFORM Solves for the polynomial coefficients for h(x)=1
%   Since the spatial shape is a uniform constant, the electric field energy
%   is exactly preserved. The energy constraint becomes purely linear.

    % 1. Create the uniform spatial shape (rank-1 tensor of all ones)
    % T_ones = ht.ones(f_ht.tree);
    T_ones = f_ht;

    % 2. Calculate global moments of the uniform basis T_ones
    M_00 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 0, 0);
    M_10 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 1, 0);
    M_01 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 0, 1);
    M_20 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 2, 0);
    M_02 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 0, 2);
    M_11 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 1, 1);
    M_30 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 3, 0);
    M_03 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 0, 3);
    M_12 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 1, 2);
    M_21 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 2, 1);
    M_40 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 4, 0);
    M_04 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 0, 4);
    M_22 = HT_compute_moments(T_ones, h_mesh, v_modes, v, 2, 2);

    % 3. Setup the EXACT FULL 4x4 linear system A*c = R
    A = zeros(4,4);
    
    % Mass row
    A(1,:) = [M_00, M_10, M_01, M_20 + M_02];
    % Momentum X row
    A(2,:) = [M_10, M_20, M_11, M_30 + M_12];
    % Momentum Y row
    A(3,:) = [M_01, M_11, M_02, M_21 + M_03];
    % Kinetic Energy row (Replaces the quadratic constraint)
    A(4,:) = 0.5 * [M_20 + M_02, M_30 + M_12, M_21 + M_03, M_40 + 2*M_22 + M_04];

    % 4. Deficit vector R
    R = [target_mass - current_mass;
         target_mom_x - current_mom_x;
         target_mom_y - current_mom_y;
         target_energy - current_energy];

    % 5. Solve directly for c (c_1, c_2, c_3, c_4)
    c = A \ R;
end