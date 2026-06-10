function [c,eta_tilde_vec,Gamma1_tilde,Gamma2_tilde] = HT_Get_coeff2(f1_ht, h_mesh, v_modes,x_modes, v, target_mass, target_mom_x, target_mom_y, target_energy, current_mass, current_mom_x, current_mom_y, current_energy, phix_ht_current, base_rel_eps, tol_ratio, min_rank, max_rank, kx)
% HT_GET_COEFF Solves for the polynomial coefficients to enforce conservation

    % 1. Calculate global moments of f1
    M_00 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 0);
    M_10 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 1, 0);
    M_01 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 1);
    M_20 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 2, 0);
    M_02 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 2);
    M_11 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 1, 1);
    M_30 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 3, 0);
    M_03 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 3);
    M_12 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 1, 2);
    M_21 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 2, 1);
    M_40 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 4, 0);
    M_04 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 0, 4);
    M_22 = HT_compute_moments(f1_ht, h_mesh, v_modes, v, 2, 2);

    % --- 2. Calculate scaled eta_tilde vector using global moments ---
    % eta_tilde = eta / eta_0
    eta_tilde_vec = [1, M_10/M_00, M_01/M_00, (M_20+M_02)/M_00];

    % --- 3. Compute Gamma1_tilde and Gamma2_tilde ---
    % Get the base density of the first mode (rho_00)
    rho_ht_00 = compute_rho(f1_ht, v_modes, h_mesh(v_modes));
    
    % Solve Poisson for rho_00
    [~, phix_ht_00, ~] = get_phi_FFT(rho_ht_00, base_rel_eps, tol_ratio, min_rank, max_rank, kx);
    
    % Spatial volume element dV
    dV = prod(h_mesh(x_modes));
    
    Gamma1_tilde = 0;
    Gamma2_tilde = 0;
    
    % Integrate over spatial dimensions
    for d = 1:length(x_modes)
        dot_phi_phi00 = ht.inner_product(phix_ht_current{d}, phix_ht_00{d});
        dot_phi00_phi00 = ht.inner_product(phix_ht_00{d}, phix_ht_00{d});
        
        Gamma1_tilde = Gamma1_tilde + dV * dot_phi_phi00;
        Gamma2_tilde = Gamma2_tilde + dV * dot_phi00_phi00;
    end
    Gamma1_tilde = real(Gamma1_tilde);
    Gamma2_tilde = real(Gamma2_tilde);
    % --- 4. Setup linear system Ac = R ---
    A = zeros(3,4);
    A(1,:) = [M_00, M_10, M_01, M_20 + M_02];
    A(2,:) = [M_10, M_20, M_11, M_30 + M_12];
    A(3,:) = [M_01, M_11, M_02, M_21 + M_03];

    R = [target_mass - current_mass;
         target_mom_x - current_mom_x;
         target_mom_y - current_mom_y];

    % --- 5. Find Particular Solution and Nullspace ---
    c_p = pinv(A) * R;      
    n = null(A);            
    n = n(:,1); 

    % --- 6. Total Energy Quadratic Equation ---
    delta_E = target_energy - current_energy;
    K_vec = 0.5 * [M_20 + M_02, M_30 + M_12, M_21 + M_03, M_40 + 2*M_22 + M_04];
    
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
        
        % Choose smallest correction
        if abs(lambda1) < abs(lambda2)
            lambda = lambda1;
        else
            lambda = lambda2;
        end
    end
    
    % Final coefficients
    c = c_p + lambda * n;
end