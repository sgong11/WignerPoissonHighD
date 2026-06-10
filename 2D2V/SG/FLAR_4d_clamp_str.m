function [f_ht_new] = FLAR_4d_clamp_str(f_ht,phi_ht,phix_ht,phixy_ht,x_extend,N_mesh,h_mesh,dom_max,dom_len,v_modes,n_ghost,dt,base_rel_eps,tol_ratio,min_rank,max_rank,H)
%FLAR_4D_CLAMP_STR Advance the Wigner force substep in 2D2V.
%
%   The force step is evaluated in velocity Fourier space with clamped
%   interpolation of the potential and its derivatives. The sampled result is
%   reconstructed in HT format by HTACA and paired with the conjugate symmetry
%   modes required for a real-valued distribution.
%
%   Inputs:
%       f_ht    - 4D distribution f(x,v,t^n) in HT format (HTD of 4D 
%                   tensor on the 2D2V grid)
%       phi_ht  - electrostatic potential in HT format
%       phix_ht - first spatial derivatives of phi in HT format
%       phixy_ht - mixed derivative information used by interpolation
%       x_extend - 1 x 2 cell arrays containing extended grid with ghost
%                   cells in x, y directions (see
%                   set_spatial_grid_info_4d.m)
%       N_mesh  - grid size in all four phase-space dimensions
%       h_mesh  - vector of grid spacing in all 4 dimensions
%       dom_max - vector of domain upper bounds in all 4 dimensions
%       dom_len - vector of domain lengths in all 4 dimensions
%       v_modes - indices of velocity dimensions
%       n_ghost - number of ghost cells used in local interpolation
%       dt      - current time step size
%       base_rel_eps - base relative tolerance for HTACA truncation
%       tol_ratio    - additional ratio factor for controlling adaptive 
%                      truncation (see ht.HTACA)
%       min_rank,max_rank       - minimum and maximum allowed ranks for 4d 
%                                 HTACA
%       H       - scaled Planck constant used in the phase factor
%
%   Outputs:
%       f_ht_new    - updated 4D distribution f(x,v,t^{n+1}) in HT format (HTD 
%                       of 4D tensor on the 2D2V grid)
    
    Nv = N_mesh(v_modes(1));
    opp_index_vx = mod(Nv-(1:Nv)+1,Nv)+1;
    vmax = dom_max(v_modes(1));
    k_mod_vx = pi * ((0:Nv-1) - Nv/2) / vmax;
    f_ht.U{5} = fftshift(fft(f_ht.U{5},f_ht.tree.modesizes(v_modes(1)),1),1);
    
    Nv = N_mesh(v_modes(2));
    opp_index_vy = mod(Nv-(1:Nv)+1,Nv)+1;
    vmax = dom_max(v_modes(2));
    k_mod_vy = pi * ((0:Nv-1) - Nv/2) / vmax;
    f_ht.U{7} = fftshift(fft(f_ht.U{7},f_ht.tree.modesizes(v_modes(2)),1),1);
    
    mid = Nv/2;
%% Keep the zero-frequency velocity modes.
    f_ht_new0   = f_ht;
    f_ht_new0.U{5}(1:mid,:) = 0;
    f_ht_new0.U{5}(mid+2:end,:) = 0;
    f_ht_new0.U{7}(1:mid,:) = 0;
    f_ht_new0.U{7}(mid+2:end,:) = 0;

 %% Reconstruct one nonzero velocity-frequency block.
    tree0 = f_ht.tree;
    tree0.modesizes(v_modes(1)) = tree0.modesizes(v_modes(1))/2 - 1;
    tree0.modesizes(v_modes(2)) = tree0.modesizes(v_modes(2))/2 + 1;
    tree0.freedom(5) = tree0.modesizes(v_modes(1));
    tree0.freedom(7) = tree0.modesizes(v_modes(2));
    tree0.freedom(2) = tree0.freedom(5)*tree0.freedom(4);
    tree0.freedom(3) = tree0.freedom(6)*tree0.freedom(7);
    tree0.freedom(1) = tree0.freedom(2)*tree0.freedom(3);
    
    fnc = @(i,j,p,q) FL_HT2D2V_handle_clamp(i,j+1,p,q,f_ht,phi_ht,phix_ht,phixy_ht,x_extend,k_mod_vx,k_mod_vy,h_mesh,N_mesh,dom_len,n_ghost,dt,H);
    max_rank = tree0.freedom; max_rank(1) = 1;
    f_ht_new = ht.HTACA(fnc,tree0,base_rel_eps,tol_ratio,min_rank,max_rank);
    
    f_ht_new.U{5}(2:mid,:) = f_ht_new.U{5}(1:mid-1,:);
    f_ht_new.U{5}(1,:) = 0;
    f_ht_new.U{5}(mid+1:N_mesh(v_modes(1)),:) = 0;
    f_ht_new.U{7}(mid+2:N_mesh(v_modes(2)),:) = 0;
    f_ht_new.tree = f_ht.tree;

    % its symmetry
    f_ht_new1 = f_ht_new;
    f_ht_new1.B{1} = conj(f_ht_new1.B{1});
    f_ht_new1.B{2} = conj(f_ht_new1.B{2});
    f_ht_new1.B{3} = conj(f_ht_new1.B{3});
    f_ht_new1.U{4} = conj(f_ht_new1.U{4});
    f_ht_new1.U{6} = conj(f_ht_new1.U{6});
    f_ht_new1.U{5} = conj(f_ht_new1.U{5}(opp_index_vx(1:N_mesh(v_modes(1))),:));
    f_ht_new1.U{7} = conj(f_ht_new1.U{7}(opp_index_vy(1:N_mesh(v_modes(2))),:));
    f_ht_new1.tree = f_ht.tree;
    f_ht_new2 = ht.add(f_ht_new,f_ht_new1);

     %% Reconstruct the complementary velocity-frequency block.
    tree0 = f_ht.tree;
    tree0.modesizes(v_modes(1)) = tree0.modesizes(v_modes(1))/2 + 1;
    tree0.modesizes(v_modes(2)) = tree0.modesizes(v_modes(2))/2 - 1;
    tree0.freedom(5) = tree0.modesizes(v_modes(1));
    tree0.freedom(7) = tree0.modesizes(v_modes(2));
    tree0.freedom(2) = tree0.freedom(5)*tree0.freedom(4);
    tree0.freedom(3) = tree0.freedom(6)*tree0.freedom(7);
    tree0.freedom(1) = tree0.freedom(2)*tree0.freedom(3);
    
    fnc = @(i,j,p,q) FL_HT2D2V_handle_clamp(i,j,p,q+mid+1,f_ht,phi_ht,phix_ht,phixy_ht,x_extend,k_mod_vx,k_mod_vy,h_mesh,N_mesh,dom_len,n_ghost,dt,H);
    max_rank = tree0.freedom; max_rank(1) = 1;
    f_ht_new = ht.HTACA(fnc,tree0,base_rel_eps,tol_ratio,min_rank,max_rank);
    
    f_ht_new.U{7}(mid+2:N_mesh(v_modes(2)),:) = f_ht_new.U{7}(1:mid-1,:);
    f_ht_new.U{5}(mid+2:N_mesh(v_modes(1)),:) = 0;
    f_ht_new.U{7}(1:mid+1,:) = 0;
    f_ht_new.tree = f_ht.tree;

    % its symmetry
    f_ht_new1 = f_ht_new;
    f_ht_new1.B{1} = conj(f_ht_new1.B{1});
    f_ht_new1.B{2} = conj(f_ht_new1.B{2});
    f_ht_new1.B{3} = conj(f_ht_new1.B{3});
    f_ht_new1.U{4} = conj(f_ht_new1.U{4});
    f_ht_new1.U{6} = conj(f_ht_new1.U{6});
    f_ht_new1.U{5} = conj(f_ht_new1.U{5}(opp_index_vx(1:N_mesh(v_modes(1))),:));
    f_ht_new1.U{7} = conj(f_ht_new1.U{7}(opp_index_vy(1:N_mesh(v_modes(2))),:));
    % f_ht_new1.U{5}(1:mid + 1,:) = 0;
    % f_ht_new1.U{7}(2:mid,:) = 0;
    f_ht_new1.tree = f_ht.tree;
    f_ht_new1 = ht.add(f_ht_new,f_ht_new1);

    f_ht_new = ht.add(f_ht_new1,f_ht_new2);
    f_ht_new = ht.add(f_ht_new,f_ht_new0);
    % [subtree_pool, ~, ~, ~, ~] = ht.get_subtree_pool(f_ht.tree);
    % [~,tol_rounding,~,~] = initialize_tolerance(subtree_pool,base_rel_eps,tol_ratio);
    % 
    % f_ht_new = ht.truncate(f_ht_new,tol_rounding(1),min_rank,max_rank); 

    f_ht_new.U{5} = ifft(ifftshift(f_ht_new.U{5},1),f_ht_new.tree.modesizes(v_modes(1)),1);
    f_ht_new.U{7} = ifft(ifftshift(f_ht_new.U{7},1),f_ht_new.tree.modesizes(v_modes(2)),1);

end
