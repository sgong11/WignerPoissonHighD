function [f_ht,E_field,E_energy,mass,momentum,energy] = initial_condition_4d(x,v,kx,dom_len,h_mesh,x_modes,v_modes,tree_4d,min_rank,max_rank,min_rank_2d,max_rank_2d,prob)
%INITIAL_CONDITION_4D Build a 2D2V initial state and its invariants.
%
%   [f_ht, E_field, E_energy, mass, momentum, energy] =
%       initial_condition_4d(x, v, kx, dom_len, h_mesh,
%                            x_modes, v_modes, tree_4d,
%                            min_rank, max_rank, min_rank_2d, max_rank_2d, prob)
%
%   This function generates the initial distribution f(x,v,0) in the
%   four-dimensional phase space (x, vx, y, vy). It solves Poisson's
%   equation for the initial electric field, normalizes the mass, and
%   returns mass, momentum, and total energy diagnostics.
%
%   Inputs:
%       x           - 1 x 2 cell array containing 1D spatial grid arrays in
%                     x and y
%       v           - 1 x 2 cell array containing 1D velocity grid arrays in
%                     vx and vy
%       kx          - 1 x 2 cell array containing fundamental frequency 
%                     vectors for FFT in x and y
%       dom_len     - domain length vector in each dimension
%       h_mesh      - grid spacing vector in each dimension
%       x_modes     - vector of mode indices corresponding to spatial dimensions
%       v_modes     - vector of mode indices corresponding to velocity dimensions 
%       tree_4d     - 4D dimension tree structure for HTACA
%       min_rank    - minimum rank allowed at each node of the 4D tree
%       max_rank    - maximum rank allowed at each node of the 4D tree
%       min_rank_2d - minimum rank allowed at each node of the 2D tree
%       max_rank_2d - maximum rank allowed at each node of the 2D tree
%       prob        - problem type identifier:
%                     1 = weak Landau damping,
%                     2 = strong Landau damping,
%                     3 = two-stream instability,
%                     4 = two-stream benchmark initialized from 1D1V,
%                     5-7 = bump-on-tail variants
%
%   Outputs:
%       f_ht     - initial distribution f(x,v,0) in HT format (HTD of 4D 
%                   tensor on the 2D2V grid)
%       E_field  - 1 x 2 cell arrays containing
%                   electric field components in HT format (HTD 
%                       of 2D tensors on the spatial grid)
%       E_energy - total initial electric field energy
%       mass     - total initial mass
%       momentum - 1 x 2 array containing initial momentum components
%       energy   - total initial energy (electric field + kinetic)

    % Set truncation parameters (special setting for initial condition)
    base_rel_eps = 1e-11;
    tol_ratio    = 0.1;
    
    if prob == 1 % Weak Landau damping
        a    = 0.01;
        data = @(i,j,p,q) initial_dis_landau_damping(i,j,p,q,x,v,x_modes,v_modes,a);
        f_ht = ht.HTACA(data,tree_4d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 2 % Strong Landau damping
        a    = 0.5;
        data = @(i,j,p,q) initial_dis_landau_damping(i,j,p,q,x,v,x_modes,v_modes,a);
        f_ht = ht.HTACA(data,tree_4d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 3 % Two-stream instability
        a   = 1e-3;
        v0  = 2.4;
        data = @(i,j,p,q) initial_dis_tsi(i,j,p,q,x,v,x_modes,v_modes,a,v0);
        f_ht = ht.HTACA(data,tree_4d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 4 % Two-stream benchmark initialized from the 1D1V profile
        data = @(i,j,p,q) initial_dis_tsi_2(i,j,p,q,x,v,x_modes,v_modes);
        f_ht = ht.HTACA(data,tree_4d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 5 % Bump-on-tail
        data = @(i,j,p,q) initial_dis_bot(i,j,p,q,x,v,x_modes,v_modes);
        f_ht = ht.HTACA(data,tree_4d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 6 % Bump-on-tail with diagonal beam
        data = @(i,j,p,q) initial_dis_bot2(i,j,p,q,x,v,x_modes,v_modes);
        f_ht = ht.HTACA(data,tree_4d,base_rel_eps,tol_ratio,min_rank,max_rank);
    elseif prob == 7 % Rotated anisotropic bump-on-tail
        base_rel_eps = 1e-7;
        data = @(i,j,p,q) initial_dis_bot3(i,j,p,q,x,v,x_modes,v_modes);
        f_ht = ht.HTACA(data,tree_4d,base_rel_eps,tol_ratio,min_rank,max_rank);
    else
        error('Unsupported problem selector: %d', prob);
    end

    % Compute the charge density of the initial distribution and get the
    % electric field via the low-rank FFT method.
    rho_ht  = compute_rho(f_ht,v_modes,h_mesh(v_modes));
    E_field = get_E_field_FFT(rho_ht,dom_len(x_modes),h_mesh(x_modes),base_rel_eps,tol_ratio,min_rank_2d,max_rank_2d,kx);
    
    % Compute electric energy, total mass, momentum components, and total energy.
    E_energy = real(compute_electric_field_energy(h_mesh(x_modes),E_field));
    mass     = compute_tmass_2D2V(f_ht,h_mesh);
    f_ht.B{1} = f_ht.B{1}*dom_len(x_modes(1))*dom_len(x_modes(2))/mass;
    mass     = compute_tmass_2D2V(f_ht,h_mesh);
    momentum = compute_momentum(f_ht,h_mesh,x_modes,v_modes,v);
    energy   = E_energy + compute_kinetic_energy(f_ht,h_mesh,x_modes,v_modes,v);

end

function f_0 = initial_dis_landau_damping(i,j,p,q,x,v,x_modes,v_modes,a)
%INITIAL_DIS_LANDAU_DAMPING Landau damping initial condition.
%
%   f_0 = initial_dis_landau_damping(i,j,p,q,x,v,x_modes,v_modes,a)
%
%   The output is vectorized for HTACA sample batches.

    coor = {i j p q};

    f_0 = 0.5/pi*( 1 + a*(cos(0.5*x{1}(coor{x_modes(1)})) + cos(0.5*x{2}(coor{x_modes(2)})) ) ) ...
                        .*exp(-0.5*v{1}(coor{v_modes(1)}).^2).*exp(-0.5*v{2}(coor{v_modes(2)}).^2);

end

function f_0 = initial_dis_tsi(i,j,p,q,x,v,x_modes,v_modes,a,v0)
%INITIAL_DIS_TSI Two-stream instability initial condition.
%
%   f_0 = initial_dis_tsi(i,j,p,q,x,v,x_modes,v_modes,a,v0)
%
%   The output is vectorized for HTACA sample batches.

    coor = {i j p q};

    v1_dis = exp(-0.5*(v{1}(coor{v_modes(1)})+v0).^2) + exp(-0.5*(v{1}(coor{v_modes(1)})-v0).^2);
    v2_dis = exp(-0.5*(v{2}(coor{v_modes(2)})+v0).^2) + exp(-0.5*(v{2}(coor{v_modes(2)})-v0).^2);

    f_0 = 0.125/pi*( 1 + a*(cos(0.2*x{1}(coor{x_modes(1)})) + cos(0.2*x{2}(coor{x_modes(2)})) ) ).*v1_dis.*v2_dis;

end

function f_0 = initial_dis_tsi_2(i,j,p,q,x,v,x_modes,v_modes)
%INITIAL_DIS_TSI_2 Two-stream benchmark initialized from a 1D1V profile.
%
%   f_0 = initial_dis_tsi_2(i,j,p,q,x,v,x_modes,v_modes)
%
%   The output is vectorized for HTACA sample batches.
    coor = {i j p q};
    
    v1_dis = (v{1}(coor{v_modes(1)})).^2.*exp(-0.5*(v{1}(coor{v_modes(1)})).^2);
    v2_dis = (v{2}(coor{v_modes(2)})).^2.*exp(-0.5*(v{2}(coor{v_modes(2)})).^2);
    f_0 = 1/(4*pi)*( 2 + (cos(0.5*x{1}(coor{x_modes(1)}))) + (cos(0.5*x{2}(coor{x_modes(2)})))).*v1_dis.*v2_dis;

end

function f_0 = initial_dis_bot(i,j,p,q,x,v,x_modes,v_modes)
%INITIAL_DIS_BOT Bump-on-tail initial condition.
%
%   f_0 = initial_dis_bot(i,j,p,q,x,v,x_modes,v_modes)
%
%   The output is vectorized for HTACA sample batches.
    coor = {i j p q};
    
    v1_dis1 = exp(-0.5*(v{1}(coor{v_modes(1)})).^2);
    v2_dis1 = exp(-0.5*(v{2}(coor{v_modes(2)})).^2);
    v1_dis2 = exp(-4*(v{1}(coor{v_modes(1)}) - 4.5).^2);
    v2_dis2 = exp(-4*(v{2}(coor{v_modes(2)})).^2);
    f_0 = 1/(2*pi)/0.925*( 1 - 0.04*(cos(0.3*x{1}(coor{x_modes(1)})) + (cos(0.3*x{2}(coor{x_modes(2)}))))).*(0.9*v1_dis1.*v2_dis1 + 0.2*v1_dis2.*v2_dis2);
    
end

function f_0 = initial_dis_bot2(i,j,p,q,x,v,x_modes,v_modes)
%INITIAL_DIS_BOT2 Bump-on-tail initial condition with a diagonal beam.
%
%   f_0 = initial_dis_bot2(i,j,p,q,x,v,x_modes,v_modes)
%
%   The output is vectorized for HTACA sample batches.
    coor = {i j p q};
 
    % Set the effective beam location
    vd = 4.5 / sqrt(2); 
    
    % Bulk plasma
    v1_dis1 = exp(-0.5*(v{1}(coor{v_modes(1)})).^2);
    v2_dis1 = exp(-0.5*(v{2}(coor{v_modes(2)})).^2);
    
    % Beam plasma (shifted to the diagonal)
    v1_dis2 = exp(-4*(v{1}(coor{v_modes(1)}) - vd).^2);
    v2_dis2 = exp(-4*(v{2}(coor{v_modes(2)}) - vd).^2);
    
    % Spatial perturbation (diagonal plane wave)
    k_sp = 0.3 / sqrt(2);
    spatial_pert = cos(k_sp*x{1}(coor{x_modes(1)}) + k_sp*x{2}(coor{x_modes(2)}));

    % Total distribution
    f_0 = 1/(2*pi)/0.925 * (1 - 0.04*spatial_pert) .* (0.9*v1_dis1.*v2_dis1 + 0.2*v1_dis2.*v2_dis2);
    
end

function f_0 = initial_dis_bot3(i,j,p,q,x,v,x_modes,v_modes)
%INITIAL_DIS_BOT3 Rotated anisotropic bump-on-tail initial condition.
%
%   f_0 = initial_dis_bot3(i,j,p,q,x,v,x_modes,v_modes)
%
%   The output is vectorized for HTACA sample batches.

    coor = {i j p q};

    % 1. Extract spatial and velocity coordinate arrays
    x1 = x{1}(coor{x_modes(1)});
    x2 = x{2}(coor{x_modes(2)});
    v1 = v{1}(coor{v_modes(1)});
    v2 = v{2}(coor{v_modes(2)});
    
    % 2. Define rotated velocity coordinates (parallel and perpendicular)
    v_par  = (v1 + v2) / sqrt(2);
    v_perp = (v1 - v2) / sqrt(2);
    
    % 3. Model parameters
    alpha   = 0.06; % Spatial perturbation amplitude
    u_b     = 4.5;    % Bump velocity (diagonal direction)
    v_tc    = 1;    % Core thermal velocity
    v_tb    = 0.5;  % Bump thermal velocity (parallel)
    v_tperp = 1;    % Bump thermal velocity (perpendicular)
    
    % 4. Bulk plasma (Core Maxwellian, Fc)
    % Fc = (1 / 2*pi) * exp(-(v_par^2 + v_perp^2) / 2)
    Fc = (1 / (2*pi)) * exp( -(v_par.^2 + v_perp.^2) / (2 * v_tc^2) );
    
    % 5. Beam plasma (Bump Maxwellian, Fb)
    % The bump is anisotropic, so we calculate it using v_par and v_perp
    % Fb = (1 / pi) * exp(-(v_par - 4)^2 / 2*(0.5)^2 - (v_perp^2) / 2)
    Fb = (1 / pi) * exp( -((v_par - u_b).^2) / (2 * v_tb^2) - (v_perp.^2) / (2 * v_tperp^2) );
    
    % 6. Spatial perturbation (diagonal plane wave)
    k_sp = 0.3 / sqrt(2); 
    spatial_pert = cos(k_sp * x1 + k_sp * x2);
    
    % 7. Total distribution f_0
    % f_0 = [1 + 1e-3 * cos(...)][0.9 * Fc + 0.1 * Fb]
    f_0 = (1 - alpha * spatial_pert) .* (0.9 * Fc + 0.1 * Fb);
end

