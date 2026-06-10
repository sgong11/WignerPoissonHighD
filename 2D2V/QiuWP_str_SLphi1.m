% === Driver for the 2D2V Wigner-Poisson solver with adaptive HT ranks ===
%
%   This script solves the Wigner-Poisson equation in four-dimensional
%   phase space (x, vx, y, vy) using Strang splitting, semi-Lagrangian
%   adaptive-rank (SLAR) transport, HTACA tensor compression, and the
%   spatially uniform conservation correction.
%   The configuration includes:
%   - Problem selection (Landau damping, two-stream, or bump-on-tail)
%   - Test type (single simulation or mesh-scaling experiment)
%   - Grid setup, tree structure, truncation settings
%   - Time stepping to final time T and saving the run data.

% restoredefaultpath
% rehash toolboxcache
% === Add necessary folders ===
addpath('Methed Related');
addpath('ht_toolbox');
addpath('SetUp Initialization');
addpath('SG');

clear all; % Clear all variables from the workspace 

% === Problem and test type selector ===
prob      = 4;  % 1 = weak Landau, 2 = strong Landau, 3-4 = two-stream, 5-6 = bump-on-tail
test_type = 1;  % 1 = single time evolution; 2 = mesh-scaling experiment

% === Problem acronym (used in file saving) ===
if prob == 1
    name_ = 'WLD';
elseif prob == 2
    name_ = 'SLD';
elseif prob == 3
    name_ = 'TSI';
elseif prob == 4
    name_ = 'TSI_test';
elseif prob == 5
    name_ = 'BOT';
elseif prob == 6
    name_ = 'BOT2';    
else
    error('Unsupported problem selector: %d', prob);
end

if test_type == 1 % single time evolution
    H = 1;
    % === Simulation settings ===
    T       = 50.;                      % Final time
    CFL     = 50.;                       % CFL number
    n_ghost = ceil(CFL) + 10;  % Ghost cells for interpolation
    % N_c = 256*[1,2,4,8];     % Optional mesh sizes for a refinement sweep
    N_c = 512;
    for l = 1:length(N_c)  
        N_mesh  = N_c(l)*ones(1,4);    % Grid size in all 4 dimensions [Nx Nvx Ny Nvy] 
                                        % (the order can change if you modify x_modes and v_modes)

        x_modes = [1 3];                    % Vector of mode indices corresponding to spatial dimensions
        v_modes = [2 4];                    % Vector of mode indices corresponding to velocity dimensions
    
        % === Construct physical grid, FFT frequencies, domain info ===
        [dom_min,dom_max,dom_len,h_mesh,x,v,x_extend,v_extend,kx] = set_spatial_grid_info_4d(N_mesh, n_ghost, x_modes, v_modes, prob);
        
        % === Construct 4D dimension tree for f(x,v) ===
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
                    0, 0; 
                    0, 0; 
                    0, 0; 
                    0, 0];
        dim2ind  = [4, 5, 6, 7];
        % The tree structure we use (binary, full, 7 nodes):
        %                (1)
        %              /     \
        %           (2)       (3)
        %          /  \       /  \
        %       (4)  (5)     (6)  (7)
        %        |    |       |    |
        %        x    vx      y    vy
        % ht.make_tree add additional fields to the tree structure for
        % convenience
        tree_4d  = ht.make_tree(children,dim2ind,N_mesh); 
    
        % === Construct 2D tree for charge density (automatically inferred) ===
        % The 4D tree and v_modes may change, so compute_rho.m is used to
        % infer the corresponding 2D density tree.
        ht_test  = ht.ones(tree_4d);
        rho_test = compute_rho(ht_test,v_modes,h_mesh(v_modes));
        tree_2d  = rho_test.tree;
    
        % === Set HTACA truncation and rank settings ===
        if prob == 1
            base_rel_eps = 1e-5;
        elseif prob == 2
            base_rel_eps = 1e-3;
        elseif prob == 3
            base_rel_eps = 1e-5;
        elseif prob == 4
            base_rel_eps = 1e-2;
        elseif prob == 5
            base_rel_eps = 1e-3;
        elseif prob == 6
            base_rel_eps = 1e-3;     
        else
            error('Unsupported problem selector: %d', prob);
        end
        tol_ratio    = 0.1;
        
        % min_rank    = ht.get_default_min_rank(tree_4d);
        min_rank    = [1,4,4,3,3,3,3];
        max_rank    = ht.get_default_max_rank(tree_4d);
        min_rank_2d = ht.get_default_min_rank(tree_2d);
        max_rank_2d = ht.get_default_max_rank(tree_2d);
        
        
        % === Initialize solution and fields ===
        [f_ht{1},E_field,E_energy(1),mass(1),momentum(1,:),energy(1)] = initial_condition_4d(x,v,kx,dom_len,h_mesh,x_modes,v_modes,tree_4d,min_rank,max_rank,min_rank_2d,max_rank_2d,prob);
        f_ht_unC = f_ht;

        % === Time-stepping loop ===    
        t(1) = 0;
        n = 1;
        while t(n) < T
            % Compute the next time step.
            [dt, t(n+1)]  = Setdt_SG(CFL,dom_max(v_modes(1)),h_mesh,t(n),T);
    
            % Run the SLAR method for one time step
            tic % record wall time per time step
            if n == 1
                f_ht_temp = f_ht{n};
            else
                fro_norm_app  = ht.norm(f_ht{n}); 
                f_ht_temp = ht.truncate(f_ht_temp,base_rel_eps*fro_norm_app,min_rank,max_rank); 
            end
           
            [f_ht_temp] = SLAR_4d(f_ht_temp,x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt/2,base_rel_eps,tol_ratio,min_rank,max_rank);
            mass_temp = compute_tmass_2D2V(f_ht_temp,h_mesh);
            f_ht_temp.B{1} = f_ht_temp.B{1}*dom_len(x_modes(1))*dom_len(x_modes(2))/mass_temp;
            rho_ht     = compute_rho(f_ht_temp,v_modes,h_mesh(v_modes));
            [phi_ht,phixm_ht,phixy_ht] = get_phi_FFT(rho_ht,base_rel_eps, tol_ratio, min_rank, max_rank, kx);
            [f_ht_temp] = FLAR_4d_clamp_str(f_ht_temp,phi_ht,phixm_ht,phixy_ht,x_extend,N_mesh,h_mesh,dom_max,dom_len,v_modes,n_ghost,dt,base_rel_eps,tol_ratio,min_rank,max_rank,H);

            fro_norm_app  = ht.norm(f_ht_temp); 
            f_ht_temp = ht.truncate(f_ht_temp,base_rel_eps*fro_norm_app,min_rank,max_rank); 
    
            Imag_B(n+1)   = imag(compute_tmass_2D2V(f_ht_temp,h_mesh));
     
            [f_ht_temp] = SLAR_4d(f_ht_temp,x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt/2,base_rel_eps,tol_ratio,min_rank,max_rank);
            
            rho_ht     = compute_rho(f_ht_temp,v_modes,h_mesh(v_modes));
            [~,phixm_ht,~] = get_phi_FFT(rho_ht,base_rel_eps, tol_ratio, min_rank, max_rank, kx);

            elec_mod = 0;
            % Get current uncorrected macroscopic quantities
            [curr_mass, curr_mom, curr_energy] = HT_Macroscope_compute(f_ht_temp, h_mesh, x_modes, v_modes, v, phixm_ht,elec_mod);
            f_ht_unC{n+1} = f_ht_temp;
            % Apply conservation correction (Uniform h(x)=1 Method)
            [f_ht{n+1}, elec_mod] = HT_Conserve_correction_uniform(f_ht_temp, h_mesh, v_modes, x_modes, v, mass(1), momentum(1,:), energy(1), curr_mass, curr_mom, curr_energy,min_rank);
            % Re-compute final macroscopic quantities to log them.
            [mass(n+1), momentum(n+1,:), energy(n+1)] = HT_Macroscope_compute(f_ht{n+1}, h_mesh, x_modes, v_modes, v, phixm_ht,elec_mod);
            
            wall_time(n) = toc; % record wall time per time step
    
            % Display progress
            disp(['t = ' num2str(t(n+1)) ', r = ['  num2str(f_ht{n+1}.rank) '], ', 'wall time =' num2str(wall_time(n))]);
    
            % Update time level 
            n = n + 1;
        end    
    % === Finalize and save ===
        computing_time = sum(wall_time);
        save(['Qiu_SLphi_' name_ '_' num2str(N_mesh(1)) '_CFL' num2str(CFL) '_T' num2str(T) '_H' num2str(H) '.mat'])
    end
elseif test_type == 2
    H = 1;
    % === Simulation settings ===
    T       = 5.;                      % Final time
    CFL     = 50.;                       % CFL number
    n_ghost = ceil(CFL) + 10;
    N_c = 30*[2,4,8,16,32,64,128]; % Mesh sizes for refinement study
    for l = 1:length(N_c)  
        N_mesh  = N_c(l)*ones(1,4);            % Grid size in all 4 dimensions [Nx Nvx Ny Nvy] 
                                            % (the order can change if you modify x_modes and v_modes)
    
        x_modes = [1 3];                    % Vector of mode indices corresponding to spatial dimensions
        v_modes = [2 4];                    % Vector of mode indices corresponding to velocity dimensions
    
        % === Construct physical grid, FFT frequencies, domain info ===
        [dom_min,dom_max,dom_len,h_mesh,x,v,x_extend,v_extend,kx] = set_spatial_grid_info_4d(N_mesh, n_ghost, x_modes, v_modes, prob);
        
        % === Construct 4D dimension tree for f(x,v) ===
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
                    0, 0; 
                    0, 0; 
                    0, 0; 
                    0, 0];
        dim2ind  = [4, 5, 6, 7];
        % The tree structure we use (binary, full, 7 nodes):
        %                (1)
        %              /     \
        %           (2)       (3)
        %          /  \       /  \
        %       (4)  (5)     (6)  (7)
        %        |    |       |    |
        %        x    vx      y    vy
        % ht.make_tree add additional fields to the tree structure for
        % convenience
        tree_4d  = ht.make_tree(children,dim2ind,N_mesh); 
    
        % === Construct 2D tree for charge density (automatically inferred) ===
        % The 4D tree and v_modes may change, so compute_rho.m is used to
        % infer the corresponding 2D density tree.
        ht_test  = ht.ones(tree_4d);
        rho_test = compute_rho(ht_test,v_modes,h_mesh(v_modes));
        tree_2d  = rho_test.tree;
    
        % === Set HTACA truncation and rank settings ===
        if prob == 1
            base_rel_eps = 1e-5;
        elseif prob == 2
            base_rel_eps = 1e-3;
        elseif prob == 3
            base_rel_eps = 1e-5;
        elseif prob == 4
            base_rel_eps = 1e-2;
        elseif prob == 5
            base_rel_eps = 1e-3;
        elseif prob == 6
            base_rel_eps = 1e-3;
        else
            error('Unsupported problem selector: %d', prob);
        end
        tol_ratio    = 0.1;
        
        min_rank    = ht.get_default_min_rank(tree_4d);
        % min_rank    = [1,4,4,3,3,3,3];
        max_rank    = ht.get_default_max_rank(tree_4d);
        % max_rank    = [1,40,40,80,80,80,80];
        min_rank_2d = ht.get_default_min_rank(tree_2d);
        max_rank_2d = ht.get_default_max_rank(tree_2d);
        % max_rank_2d = [1,50,50];
        
        
        % === Initialize solution and fields ===
        [f_ht{1},E_field,E_energy(1),mass(1),momentum(1,:),energy(1)] = initial_condition_4d(x,v,kx,dom_len,h_mesh,x_modes,v_modes,tree_4d,min_rank,max_rank,min_rank_2d,max_rank_2d,prob);
    
        % === Time-stepping loop ===    
        t(1) = 0;
        n = 1;
        while t(n) < T
            % Compute the next time step.
            % [dt, t(n+1)]  = Setdt_SG(CFL,dom_max(v_modes(1)),h_mesh,t(n),T);
            dt = 0.1;
            t(n+1) = t(n)+dt;
            if t(n+1) > T
                dt = T - t(n);
                t(n+1) = T;
            end
            % Run the SLAR method for one time step
            tic % record wall time per time step
            [f_ht_temp] = SLAR_4d(f_ht{n},x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt/2,base_rel_eps,tol_ratio,min_rank,max_rank);
            rho_ht     = compute_rho(f_ht_temp,v_modes,h_mesh(v_modes));
            [phi_ht,phix_ht,phixy_ht] = get_phi_FFT(rho_ht,base_rel_eps, tol_ratio, min_rank, max_rank, kx);
            [f_ht_temp] = FLAR_4d_clamp_str(f_ht_temp,phi_ht,phix_ht,phixy_ht,x_extend,N_mesh,h_mesh,dom_max,dom_len,v_modes,n_ghost,dt,base_rel_eps,tol_ratio,min_rank,max_rank,H);
    
            Imag_B(n+1)   = imag(compute_tmass_2D2V(f_ht_temp,h_mesh));
     
            [f_ht{n+1}] = SLAR_4d(f_ht_temp,x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt/2,base_rel_eps,tol_ratio,min_rank,max_rank);
            mass(n+1)     = compute_tmass_2D2V(f_ht{n+1},h_mesh);
            f_ht{n+1}.B{1} = f_ht{n+1}.B{1}*dom_len(x_modes(1))*dom_len(x_modes(2))/mass(n+1);
            mass(n+1)     = compute_tmass_2D2V(f_ht{n+1},h_mesh);
            
            wall_time(n) = toc; % record wall time per time step
    
            % Display progress
            % disp(['t = ' num2str(t(n+1)) ', r = ['  num2str(f_ht{n+1}.rank) '], ', 'wall time =' num2str(wall_time(n))]);
    
            % Update time level 
            n = n + 1;
        end 
        computing_time(l) = sum(wall_time);
    end
    % === Finalize and save ===
    % computing_time = sum(wall_time);
    save(['Test_' name_ '_' num2str(N_mesh(1)) '_CFL' num2str(CFL) '_T' num2str(T) '_H' num2str(H) '.mat'])
else
    error('Unsupported test_type: %d', test_type);
end
