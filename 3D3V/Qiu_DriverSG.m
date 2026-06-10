% === Main driver for the 3D3V Wigner-Poisson adaptive-rank solver ===
%
%   This script runs the 3D3V Wigner-Poisson equation with semi-Lagrangian
%   adaptive-rank updates, Hierarchical Tucker tensor compression, and the
%   Fourier-space velocity update used in the manuscript experiments.

% === Add necessary folders ===
addpath('Methed Related');
addpath('ht_toolbox');
addpath('SetUp Initialization');
addpath('Gemini');

clear all; % Clear all variables from the workspace 

% === Problem selector ===
prob      = 4;  % 1 = Weak Landau damping; 2 = Strong Landau damping; 3 = two-stream instability; 4 = Wigner-Poisson two-stream case
test_type = 1;  % Retained driver path: single time evolution

% === Problem acronym (used in file saving) ===
if prob == 1
    name_ = 'WLD';
elseif prob == 2
    name_ = 'SLD';
elseif prob == 3
    name_ = 'TSI';
elseif prob == 4
    name_ = 'TSI_WP';
end

if test_type == 1
    H = 8;
    % === Simulation settings ===
    T       = 30.;                      % Final time
    CFL     = 5.;                       % CFL number
    n_ghost = ceil(CFL) + 10;           % Number of ghost cells for interpolation
    N_mesh  = 256*ones(1,6);            % Grid size in all 6 dimensions [Nx Nvx Ny Nvy Nz Nvz] 
                                        % (the order can change if you modify x_modes and v_modes)

    x_modes = [1 3 5];                  % Vector of mode indices corresponding to spatial dimensions
    v_modes = [2 4 6];                  % Vector of mode indices corresponding to velocity dimensions

    % === Construct physical grid, FFT frequencies, and domain info ===
    [dom_min,dom_max,dom_len,h_mesh,x,v,x_extend,v_extend,kx] = set_spatial_grid_info_6d(N_mesh, n_ghost, x_modes, v_modes, prob);
    
    % === Construct 6D dimension tree for f(x,v) ===
    % Must use ht.make_tree to generate a valid full tree structure.
    %
    % The tree is defined by the 'children' array:
    %   - children(i,:) = [l, r] means node i has left child l and right child r.
    %   - children(i,:) = [0, 0] means node i is a leaf node (i.e., corresponds to a physical dimension).
    %
    % The physical dimensions are mapped to leaf node indices by 'dim2ind'
    %
    % We enforce the root node index to be 1 in order to use the ht_toolbox
    children = [2, 3; 
                4, 5; 
                6, 7; 
                8, 9; 
                10, 11; 
                0, 0; 
                0, 0; 
                0, 0; 
                0, 0; 
                0, 0; 
                0, 0];
    dim2ind  = [8, 9, 10, 11, 6, 7];
    % The tree structure we use (binary, full, 11 nodes):
    %                (1)
    %              /     \
    %           (2)       (3)
    %          /  \       /  \
    %       (4)  (5)     (6)  (7)
    %      / \   / \      |    | 
    %    (8)(9) (10)(11)  |    |        
    %     |  |   |   |    |    |
    %     x vx   y  vy    z    vz
    % ht.make_tree add additional fields to the tree structure for
    % convenience
    tree_6d  = ht.make_tree(children,dim2ind,N_mesh); 

    % === Construct 3D tree for charge density (automatically inferred) ===
    % The 6D tree and v_modes may change, so compute_rho infers the
    % corresponding 3D tree automatically.
    ht_test  = ht.ones(tree_6d);
    rho_test = compute_rho(ht_test,v_modes,h_mesh(v_modes));
    tree_3d  = rho_test.tree;

    % === Set HTACA truncation and rank settings ===
    if prob == 1
        base_rel_eps = 1e-5;
    elseif prob == 2
        base_rel_eps = 1e-4;
    elseif prob == 3
        base_rel_eps = 1e-3;
    elseif prob == 4
        base_rel_eps = 1e-2;
    end
    tol_ratio    = 0.1;
    if prob ~= 4
        min_rank    = ht.get_default_min_rank(tree_6d);
        max_rank    = ht.get_default_max_rank(tree_6d);
        min_rank_3d = ht.get_default_min_rank(tree_3d);
        max_rank_3d = ht.get_default_max_rank(tree_3d);
    else
        min_rank    = [1 1 1 1 1 3 3 3 3 3 3];
        max_rank    = ht.get_default_max_rank(tree_6d);
        min_rank_3d = [1 1 3 3 3];
        max_rank_3d = ht.get_default_max_rank(tree_3d);
    end

    % === Initialize solution and fields ===
    [f_ht{1},E_field,E_energy(1),mass(1),momentum(1,:),energy(1)] = initial_condition_6d2(x,v,kx,dom_len,h_mesh,x_modes,v_modes,tree_6d,min_rank,max_rank,min_rank_3d,max_rank_3d,prob);

    % === Time-stepping loop ===    
    t(1) = 0;
    n = 1;
    while t(n) < T
        % Compute the next time step and time level.
        [dt, t(n+1)]  = Setdt_SG(CFL,dom_max(v_modes(1)),h_mesh,t(n),T);

        % Run the splitting update for one time step.
        tic % record wall time per time step
        if n == 1
            f_ht_temp = f_ht{n};
        else
            fro_norm_app  = ht.norm(f_ht{n}); 
            f_ht_temp = ht.truncate(f_ht_temp,base_rel_eps*fro_norm_app,min_rank,max_rank); 
        end
        
        [f_ht_temp] = AI_SLAR_6d_WP(f_ht_temp,x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt/2,base_rel_eps,tol_ratio,min_rank,max_rank);
        mass_temp = compute_tmass(f_ht_temp,h_mesh);
        f_ht_temp.B{1} = f_ht_temp.B{1}*dom_len(x_modes(1))*dom_len(x_modes(2))*dom_len(x_modes(3))/mass_temp;
        rho_ht     = compute_rho(f_ht_temp,v_modes,h_mesh(v_modes));
        
        % Compute the potential, gradient, and mixed derivatives.
        [phi_ht, phix_ht, phi_cross] = get_phi_FFT(rho_ht,base_rel_eps, tol_ratio, min_rank, max_rank, kx);
        
        % Apply the Fourier-space velocity update.
        [f_ht_temp] = AI_FLAR_6d_str(f_ht_temp,phi_ht,phix_ht,phi_cross,x_extend,N_mesh,h_mesh,dom_max,x_modes,v_modes,n_ghost,dt,base_rel_eps,tol_ratio,min_rank,max_rank,H);
        
        Imag_B(n+1)   = imag(compute_tmass(f_ht_temp,h_mesh));

        fro_norm_app  = ht.norm(f_ht_temp); 
        f_ht_temp = ht.truncate(f_ht_temp,base_rel_eps*fro_norm_app,min_rank,max_rank); 

        [f_ht_temp] = AI_SLAR_6d_WP(f_ht_temp,x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt/2,base_rel_eps,tol_ratio,min_rank,max_rank);
        rho_ht     = compute_rho(f_ht_temp,v_modes,h_mesh(v_modes));
        
        % Only the gradient is needed for macroscopic quantities and correction.
        [~, phixm_ht, ~] = get_phi_FFT(rho_ht,base_rel_eps, tol_ratio, min_rank, max_rank, kx);
        
        
        elec_mod = 0;
        % Get current uncorrected macroscopic quantities
        [curr_mass, curr_mom, curr_energy] = HT_Macroscope_compute(f_ht_temp, h_mesh, x_modes, v_modes, v, phixm_ht,elec_mod);
        
        % Apply conservation correction
        [f_ht{n+1},elec_mod] = HT_Conserve_correction3D_uniform(f_ht_temp, h_mesh, v_modes, x_modes,v, mass(1), momentum(1,:), energy(1), curr_mass, curr_mom, curr_energy, phixm_ht,base_rel_eps,tol_ratio, min_rank, max_rank,kx);
        % Re-compute final macroscopic quantities to log them
        [mass(n+1), momentum(n+1,:), energy(n+1)] = HT_Macroscope_compute(f_ht{n+1}, h_mesh, x_modes, v_modes, v, phixm_ht,elec_mod);
        
        wall_time(n) = toc; % record wall time per time step

        % Display progress
        disp(['t = ' num2str(t(n+1)) ', r = ['  num2str(f_ht{n+1}.rank) '], ', 'wall time =' num2str(wall_time(n))]);

        % Update time level 
        n = n + 1;
    end    

    % === Finalize and save ===
    computing_time = sum(wall_time);
    save(['Qiu_' name_ '_' num2str(N_mesh(1)) '_CFL' num2str(CFL) '_T' num2str(T) '_H' num2str(H) '.mat'])

end
