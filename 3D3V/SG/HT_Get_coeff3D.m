function [c, eta_tilde_vec, Gamma1_tilde, Gamma2_tilde] = HT_Get_coeff3D(f1_ht, h_mesh, v_modes, x_modes, v, target_mass, target_mom, target_energy, current_mass, current_mom, current_energy, phix_ht_current, base_rel_eps, tol_ratio, min_rank, max_rank, kx)
% HT_GET_COEFF3D Solve polynomial correction coefficients for 3D3V conservation

    % 1. Calculate global velocity moments of the correction basis.
    % Mass
    M_000 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 0, 0);
    
    % Momentum
    M_100 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 1, 0, 0);
    M_010 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 1, 0);
    M_001 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 0, 1);
    
    % Kinetic-energy moments
    M_200 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 2, 0, 0);
    M_020 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 2, 0);
    M_002 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 0, 2);
    
    % Mixed momentum moments
    M_110 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 1, 1, 0);
    M_101 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 1, 0, 1);
    M_011 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 1, 1);
    
    % Third-order moments
    M_300 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 3, 0, 0);
    M_030 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 3, 0);
    M_003 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 0, 3);
    
    M_210 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 2, 1, 0);
    M_201 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 2, 0, 1);
    M_120 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 1, 2, 0);
    M_021 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 2, 1);
    M_102 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 1, 0, 2);
    M_012 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 1, 2);
    
    % Fourth-order moments
    M_400 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 4, 0, 0);
    M_040 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 4, 0);
    M_004 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 0, 4);
    
    M_220 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 2, 2, 0);
    M_202 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 2, 0, 2);
    M_022 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 2, 2);

    % --- 2. Calculate scaled eta_tilde vector using global moments ---
    % eta_tilde = eta / eta_0
    eta_tilde_vec = [1, M_100/M_000, M_010/M_000, M_001/M_000, (M_200+M_020+M_002)/M_000];

    % --- 3. Compute Gamma1_tilde and Gamma2_tilde ---
    % Get the density of the correction basis (rho_000).
    rho_ht_000 = compute_rho(f1_ht, v_modes, h_mesh(v_modes));
    
    % Solve Poisson for rho_000
    [~, phix_ht_000] = get_phi_FFT(rho_ht_000, base_rel_eps, tol_ratio, min_rank, max_rank, kx);
    
    % Spatial volume element dV
    dV = prod(h_mesh(x_modes));
    
    Gamma1_tilde = 0;
    Gamma2_tilde = 0;
    
    % Integrate over spatial dimensions
    for d = 1:length(x_modes)
        dot_phi_phi00 = ht.inner_product(phix_ht_current{d}, phix_ht_000{d});
        dot_phi00_phi00 = ht.inner_product(phix_ht_000{d}, phix_ht_000{d});
        
        Gamma1_tilde = Gamma1_tilde + dV * dot_phi_phi00;
        Gamma2_tilde = Gamma2_tilde + dV * dot_phi00_phi00;
    end
    Gamma1_tilde = real(Gamma1_tilde);
    Gamma2_tilde = real(Gamma2_tilde);

    % --- 4. Setup linear system Ac = R ---
    A = zeros(4,5);
    A(1,:) = [M_000, M_100, M_010, M_001, M_200 + M_020 + M_002];
    A(2,:) = [M_100, M_200, M_110, M_101, M_300 + M_120 + M_102];
    A(3,:) = [M_010, M_110, M_020, M_011, M_210 + M_030 + M_012];
    A(4,:) = [M_001, M_101, M_011, M_002, M_201 + M_021 + M_003];

    R = [target_mass - current_mass;
         target_mom(1) - current_mom(1);
         target_mom(2) - current_mom(2);
         target_mom(3) - current_mom(3)];

    % --- 5. Find a particular solution and nullspace direction ---
    c_p = pinv(A) * R;      
    n = null(A);            
    n = n(:,1); 

    % --- 6. Total Energy Quadratic Equation ---
    delta_E = target_energy - current_energy;
    
    K_vec = 0.5 * [M_200 + M_020 + M_002, ...
                   M_300 + M_120 + M_102, ...
                   M_210 + M_030 + M_012, ...
                   M_201 + M_021 + M_003, ...
                   M_400 + M_040 + M_004 + 2*(M_220 + M_202 + M_022)];
    
    K_p = dot(K_vec, c_p);
    K_n = dot(K_vec, n);
    
    % Use eta_tilde instead of eta
    eta_tilde_p = dot(eta_tilde_vec, c_p);
    eta_tilde_n = dot(eta_tilde_vec, n);
    
    % Formulate quadratic coefficients (a*lam^2 + b*lam + c_quad = 0)
    a = 0.5 * (eta_tilde_n^2) * Gamma2_tilde;
    b = K_n + eta_tilde_n * Gamma1_tilde + eta_tilde_p * eta_tilde_n * Gamma2_tilde;
    c_quad = K_p + eta_tilde_p * Gamma1_tilde + 0.5 * (eta_tilde_p^2) * Gamma2_tilde - delta_E;

    % --- 7. Solve for lambda ---
    if abs(a) < 1e-14
        lambda = -c_quad / b; % Frozen-field fallback
    else
        discriminant = b^2 - 4*a*c_quad;
        if discriminant < 0
            error('No real solution for energy constraint correction.');
        end
        lambda1 = (-b + sqrt(discriminant)) / (2*a);
        lambda2 = (-b - sqrt(discriminant)) / (2*a);
        
            % Choose the smallest correction.
        if abs(lambda1) < abs(lambda2)
            lambda = lambda1;
        else
            lambda = lambda2;
        end
    end
    
    % Final coefficients
    c = c_p + lambda * n;
end
