function [f_ht,E_field,E_energy,mass,momentum,energy] = initial_condition_6d2(x,v,kx,dom_len,h_mesh,x_modes,v_modes,tree_6d,min_rank,max_rank,min_rank_3d,max_rank_3d,prob)
%INITIAL_CONDITION_6D2 Set up the initial condition for the 3D3V Wigner-Poisson problem
%and compute the corresponding initial macroscopic quantities.
%
%   [f_ht, E_field, E_energy, mass, momentum, energy] =
%       initial_condition_6d2(x, v, kx, dom_len, h_mesh,
%                            x_modes, v_modes, tree_6d,
%                            min_rank, max_rank, min_rank_3d, max_rank_3d, prob)
%
%   This function generates the initial 6D phase-space distribution \( f(x,v) \)
%   and computes the electric field using Poisson's equation. It also outputs
%   the initial mass, momentum, and total energy (kinetic + electric field).
%
%   Inputs:
%       x           - 1 x 3 cell array containing 1D spatial grid arrays in
%                     x, y, and z directions
%       v           - 1 x 3 cell array containing 1D spatial grid arrays in
%                     vx, vy, and vz directions
%       kx          - 1 x 3 cell array containing fundamental frequency 
%                     vectors for FFT (in x, y, z)
%       dom_len     - domain length vector in each dimension
%       h_mesh      - grid spacing vector in each dimension
%       x_modes     - vector of mode indices corresponding to spatial dimensions
%       v_modes     - vector of mode indices corresponding to velocity dimensions 
%       tree_6d     - 6D dimension tree structure for 6D HTACA
%       min_rank    - minimum rank allowed at each node of the 6D tree
%       max_rank    - maximum rank allowed at each node of the 6D tree
%       min_rank_3d - minimum rank allowed at each node of the 3D tree
%       max_rank_3d - maximum rank allowed at each node of the 3D tree
%       prob        - problem type identifier:
%                     1 = weak Landau damping,
%                     2 = strong Landau damping,
%                     3 = two-stream instability,
%                     4 = Wigner-Poisson two-stream case,
%                     5 = bump-on-tail instability
%
%   Outputs:
%       f_ht     - initial distribution f(x,v,0) in HT format (HTD of 6D 
%                   tensor discretized at cell centers)
%       E_field  - 1 x 3 cell arrays containing
%                   electric field components in HT format (HTD 
%                       of 3D tensors discretized at cell centers)
%       E_energy - total initial electric field energy
%       mass     - total initial mass
%       momentum - 1 x 3 array containing initial momentum components
%       energy   - total initial energy (electric field + kinetic)

    % Set truncation parameters (special setting for initial condition)
    base_rel_eps = 1e-11;
    tol_ratio    = 0.1;
    
    if prob == 1 % Weak Landau damping
        % Set up the function handle for the weak Landau damping distribution.
        a    = 0.01;
        data = @(i,j,p,q,m,n) initial_dis_landau_damping(i,j,p,q,m,n,x,v,x_modes,v_modes,a);
        % Use ht.HTACA to form the initial distribution in HT format
        f_ht = ht.HTACA(data,tree_6d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 2 % Strong Landau damping
        a    = 1./3;
        data = @(i,j,p,q,m,n) initial_dis_landau_damping(i,j,p,q,m,n,x,v,x_modes,v_modes,a);
        f_ht = ht.HTACA(data,tree_6d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 3 % Two-stream instability
        a   = 1e-3;
        v0  = 2.4;
        data = @(i,j,p,q,m,n) initial_dis_tsi(i,j,p,q,m,n,x,v,x_modes,v_modes,a,v0);
        f_ht = ht.HTACA(data,tree_6d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 4 % Wigner-Poisson two-stream case
        data = @(i,j,p,q,m,n) initial_dis_tsi2(i,j,p,q,m,n,x,v,x_modes,v_modes);
        f_ht = ht.HTACA(data,tree_6d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 5 % Bump-on-tail instability
        data = @(i,j,p,q,m,n) initial_dis_bot(i,j,p,q,m,n,x,v,x_modes,v_modes);
        f_ht = ht.HTACA(data,tree_6d,base_rel_eps,tol_ratio,min_rank,max_rank);
    end

    % Compute the charge density of the initial distribution and get the
    % electric field via the low-rank FFT method.
    rho_ht  = compute_rho(f_ht,v_modes,h_mesh(v_modes));
    E_field = get_E_field_FFT(rho_ht,dom_len(x_modes),h_mesh(x_modes),base_rel_eps,tol_ratio,min_rank_3d,max_rank_3d,kx);
    
    % Compute electric energy, total mass, momentum components, and total
    % energy
    E_energy = real(compute_electric_field_energy(h_mesh(x_modes),E_field));
    mass     = compute_tmass(f_ht,h_mesh);
    momentum = compute_momentum(f_ht,h_mesh,x_modes,v_modes,v);
    energy   = E_energy + compute_kinetic_energy(f_ht,h_mesh,x_modes,v_modes,v);

end

function f_0 = initial_dis_landau_damping(i,j,p,q,m,n,x,v,x_modes,v_modes,a)
%INITIAL_DIS_LANDAU_DAMPING Function handle for Landau damping initial condition
%
%   f_0 = initial_dis_landau_damping(i,j,p,q,m,n,x,v,x_modes,v_modes,a)
%
%   f_0 has to be a column vector if i,..n or one of them are vectors
%   reshape(f_0,[],1) if needed;

    coor = {i j p q m n};

    f_0 = 1/(2*pi)^(1.5)*( 1 + a*(cos(0.5*x{1}(coor{x_modes(1)})) + cos(0.5*x{2}(coor{x_modes(2)})) + cos(0.5*x{3}(coor{x_modes(3)})) ) ) ...
                        .*exp(-0.5*v{1}(coor{v_modes(1)}).^2).*exp(-0.5*v{2}(coor{v_modes(2)}).^2).*exp(-0.5*v{3}(coor{v_modes(3)}).^2);

end

function f_0 = initial_dis_tsi(i,j,p,q,m,n,x,v,x_modes,v_modes,a,v0)
%INITIAL_DIS_TSI Function handle for two-stream instability initial condition
%
%   f_0 = initial_dis_tsi(i,j,p,q,m,n,x,v,x_modes,v_modes,a,v0)
%
%   f_0 has to be a column vector if i,..n or one of them are vectors
%   reshape(f_0,[],1) if needed;

    coor = {i j p q m n};

    v1_dis = exp(-0.5*(v{1}(coor{v_modes(1)})+v0).^2) + exp(-0.5*(v{1}(coor{v_modes(1)})-v0).^2);
    v2_dis = exp(-0.5*(v{2}(coor{v_modes(2)})+v0).^2) + exp(-0.5*(v{2}(coor{v_modes(2)})-v0).^2);
    v3_dis = exp(-0.5*(v{3}(coor{v_modes(3)})+v0).^2) + exp(-0.5*(v{3}(coor{v_modes(3)})-v0).^2);

    f_0 = 0.125/(2*pi)^(1.5)*( 1 + a*(cos(0.2*x{1}(coor{x_modes(1)})) + cos(0.2*x{2}(coor{x_modes(2)})) + cos(0.2*x{3}(coor{x_modes(3)})) ) ).*v1_dis.*v2_dis.*v3_dis;

end

function f_0 = initial_dis_tsi2(i,j,p,q,m,n,x,v,x_modes,v_modes)
%INITIAL_DIS_TSI2 Function handle for the Wigner-Poisson two-stream initial condition
%
%   f_0 = initial_dis_tsi2(i,j,p,q,m,n,x,v,x_modes,v_modes)
%
%   f_0 has to be a column vector if i,..n or one of them are vectors
%   reshape(f_0,[],1) if needed;

    coor = {i j p q m n};

    % v1_dis = exp(-0.5*(v{1}(coor{v_modes(1)})+v0).^2) + exp(-0.5*(v{1}(coor{v_modes(1)})-v0).^2);
    % v2_dis = exp(-0.5*(v{2}(coor{v_modes(2)})+v0).^2) + exp(-0.5*(v{2}(coor{v_modes(2)})-v0).^2);
    % v3_dis = exp(-0.5*(v{3}(coor{v_modes(3)})+v0).^2) + exp(-0.5*(v{3}(coor{v_modes(3)})-v0).^2);
    v1_dis = (v{1}(coor{v_modes(1)})).^2.*exp(-0.5*(v{1}(coor{v_modes(1)})).^2);
    v2_dis = (v{2}(coor{v_modes(2)})).^2.*exp(-0.5*(v{2}(coor{v_modes(2)})).^2);
    v3_dis = (v{3}(coor{v_modes(3)})).^2.*exp(-0.5*(v{3}(coor{v_modes(3)})).^2);
    
    f_0 = 1/sqrt(2*pi)/(4*pi)*( 2 + (cos(0.5*x{1}(coor{x_modes(1)}))) + (cos(0.5*x{2}(coor{x_modes(2)})))+(cos(0.5*x{3}(coor{x_modes(3)})))).*v1_dis.*v2_dis.*v3_dis;
    
end

function f_0 = initial_dis_bot(i,j,p,q,m,n,x,v,x_modes,v_modes)
%INITIAL_DIS_BOT Function handle for the 3D bump-on-tail instability initial condition
%
%   f_0 = initial_dis_bot(i, j, p, q, m, n, x, v, x_modes, v_modes)
%
%   f_0 has to be a column vector if i,..n or one of them are vectors
%   reshape(f_0,[],1) if needed;

    % 6D coordinates: 3 spatial (i, p, m), 3 velocity (j, q, n)
    coor = {i, j, p, q, m, n};
    
    % Background Maxwellian (Variance = 1, Mean = 0 in all 3 directions)
    v1_dis1 = exp(-0.5*(v{1}(coor{v_modes(1)})).^2);
    v2_dis1 = exp(-0.5*(v{2}(coor{v_modes(2)})).^2);
    v3_dis1 = exp(-0.5*(v{3}(coor{v_modes(3)})).^2);
    
    % Bump Maxwellian (Narrower, drifting at v=4.5 in the v1 direction)
    v1_dis2 = exp(-4*(v{1}(coor{v_modes(1)}) - 4.5).^2);
    v2_dis2 = exp(-4*(v{2}(coor{v_modes(2)})).^2);
    v3_dis2 = exp(-4*(v{3}(coor{v_modes(3)})).^2);
    
    % 3D Spatial perturbation
    spatial_pert = 1 - 0.04*(cos(0.3*x{1}(coor{x_modes(1)})) + ...
                             cos(0.3*x{2}(coor{x_modes(2)})) + ...
                             cos(0.3*x{3}(coor{x_modes(3)})));
                         
    % Exact analytical normalization constant for 3D velocity integral
    % N = int(0.9 * exp(-0.5*v^2) + 0.2 * exp(-4*v^2)) d^3v
    norm_const = 0.9 * (2*pi)^1.5 + 0.2 * (pi/4)^1.5;
    
    % Calculate the final 6D initial distribution
    f_0 = (1 / norm_const) * spatial_pert .* ...
          (0.9 * v1_dis1 .* v2_dis1 .* v3_dis1 + ...
           0.2 * v1_dis2 .* v2_dis2 .* v3_dis2);
           
end
