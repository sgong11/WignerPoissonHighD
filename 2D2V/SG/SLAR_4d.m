function [f_ht_new] = SLAR_4d(f_ht,x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt,base_rel_eps,tol_ratio,min_rank,max_rank)
%SLAR_4D Advance the spatial-transport substep in 2D2V.
%
%   The substep solves x_t + v·grad_x f = 0 over the supplied time interval
%   by evaluating a WENO semi-Lagrangian reconstruction through HTACA. The
%   returned tensor is recompressed using the requested rank bounds.
%
%   Inputs:
%       f_ht    - 4D distribution f(x,v,t^n) in HT format (HTD of 4D 
%                   tensor on the 2D2V grid)
%       x_extend - 1 x 2 cell arrays containing extended grid with ghost
%                   cells in x, y directions (see
%                   set_spatial_grid_info_4d.m)
%       v       - 1 x 2 cell arrays containing velocity grid points in vx, vy
%       h_mesh  - vector of grid spacing in all 4 dimensions
%       dom_min - vector of domain lower bounds in all 4 dimensions
%       x_modes - indices of spatial dimensions
%       v_modes - indices of velocity dimensions
%       n_ghost - number of ghost cells used in local interpolation
%       dt      - current time step size
%       base_rel_eps - base relative tolerance for HTACA truncation
%       tol_ratio    - additional ratio factor for controlling adaptive 
%                      truncation (see ht.HTACA)
%       min_rank,max_rank       - minimum and maximum allowed ranks for 4d 
%                                 HTACA
%       min_rank_2d,max_rank_2d - minimum and maximum allowed ranks for 2d 
%                                 HTACA
%
%   Outputs:
%       f_ht_new    - updated 4D distribution f(x,v,t^{n+1}) in HT format (HTD 
%                       of 4D tensor on the 2D2V grid)

    % Build the sampled transport map used by HTACA.
    fnc = @(i,j,p,q) SL_HT2D2V_handle_vec(i,j,p,q,f_ht,x_extend,v,h_mesh,dom_min,x_modes,v_modes,n_ghost,dt,1e-6);
    % Apply HTACA to construct the intermediate solution in HT format.
    f_ht_new = ht.HTACA(fnc,f_ht.tree,base_rel_eps,tol_ratio,min_rank,max_rank);

end
