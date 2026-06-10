function [f_ht_new] = AI_SLAR_6d_WP(f_ht,x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt,base_rel_eps,tol_ratio,min_rank,max_rank)
%AI_SLAR_6D_WP Apply the adaptive-rank spatial transport update for 3D3V Wigner-Poisson.
%
%   This function advances the spatial-advection part of the split 3D3V
%   Wigner-Poisson update. The local semi-Lagrangian WENO interpolation is
%   exposed as a function handle and recompressed by HTACA.
%
%   Inputs:
%       f_ht    - 6D distribution f(x,v,t^n) in HT format (HTD of 6D 
%                   tensor discretized at cell centers)
%       x_extend - 1 x 3 cell arrays containing extended grid with ghost
%                   cells in x, y and z directions (see
%                   set_spatial_grid_info_6d.m)
%       v       - 1 x 3 cell array containing velocity grids in vx, vy and vz
%       h_mesh  - vector of grid spacing in all 6 dimensions
%       dom_min - vector of domain lower bounds; retained for interface consistency
%       x_modes - indices of spatial dimensions
%       v_modes - indices of velocity dimensions
%       n_ghost - number of ghost cells used in local interpolation
%       dt      - current time step size
%       base_rel_eps - base relative tolerance for HTACA truncation
%       tol_ratio    - additional ratio factor for controlling adaptive 
%                      truncation (see ht.HTACA)
%       min_rank,max_rank - minimum and maximum allowed ranks for 6D HTACA
%
%   Outputs:
%       f_ht_new - updated 6D distribution in HT format

    fnc = @(i,j,p,q,l,m) SL_HT3D3V_handle_vec(i,j,p,q,l,m,f_ht,x_extend,v,h_mesh,x_modes,v_modes,n_ghost,dt,1e-6);
    % Apply HTACA to construct the transported solution in HT format.
    f_ht_new = ht.HTACA(fnc,f_ht.tree,base_rel_eps,tol_ratio,min_rank,max_rank);
end
