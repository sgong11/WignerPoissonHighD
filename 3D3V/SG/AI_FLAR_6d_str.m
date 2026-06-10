function [f_ht_new] = AI_FLAR_6d_str(f_ht, phi_ht, phix_ht, phi_cross, x_extend, N_mesh, h_mesh, dom_max, x_modes, v_modes, n_ghost, dt, base_rel_eps, tol_ratio, min_rank, max_rank, H)
%AI_FLAR_6D_STR Apply the Fourier-space velocity update for the 3D3V Wigner-Poisson split step.
%
%   This function advances the velocity part of the split Wigner-Poisson
%   update. Velocity bases are transformed to Fourier space, the phase
%   induced by phi is evaluated through clamped 3D Hermite interpolation,
%   and the result is recompressed by HTACA.
%
%   Inputs:
%       f_ht      - 6D distribution f(x,v,t^n) in HT format
%       phi_ht    - electric potential in 3D HT format
%       phix_ht   - cell array of potential gradients in 3D HT format
%       phi_cross - mixed spectral derivatives used by 3D Hermite interpolation
%       x_extend  - 1 x 3 cell arrays containing extended spatial grids
%       N_mesh    - grid size in all 6 dimensions
%       h_mesh    - vector of grid spacing in all 6 dimensions
%       dom_max   - vector of domain upper bounds in all 6 dimensions
%       x_modes - indices of spatial dimensions
%       v_modes - indices of velocity dimensions
%       n_ghost - number of ghost cells used in local interpolation
%       dt      - current time step size
%       base_rel_eps - base relative tolerance for HTACA truncation
%       tol_ratio    - additional ratio factor for controlling adaptive 
%                      truncation (see ht.HTACA)
%       min_rank,max_rank - minimum and maximum allowed ranks for 6D HTACA
%       H       - scaled Planck parameter in the Wigner potential term
%
%   Outputs:
%       f_ht_new - updated 6D distribution in HT format
    
    Nv = N_mesh(v_modes(1));
    opp_index_vx = mod(Nv-(1:Nv)+1,Nv)+1;
    vmax = dom_max(v_modes(1));
    k_mod_vx = pi * ((0:Nv-1) - Nv/2) / vmax;
    f_ht.U{9} = fftshift(fft(f_ht.U{9},f_ht.tree.modesizes(v_modes(1)),1),1);
    
    Nv = N_mesh(v_modes(2));
    opp_index_vy = mod(Nv-(1:Nv)+1,Nv)+1;
    vmax = dom_max(v_modes(2));
    k_mod_vy = pi * ((0:Nv-1) - Nv/2) / vmax;
    f_ht.U{11} = fftshift(fft(f_ht.U{11},f_ht.tree.modesizes(v_modes(2)),1),1);
    
    Nv = N_mesh(v_modes(3));
    opp_index_vz = mod(Nv-(1:Nv)+1,Nv)+1;
    vmax = dom_max(v_modes(3));
    k_mod_vz = pi * ((0:Nv-1) - Nv/2) / vmax;
    f_ht.U{7} = fftshift(fft(f_ht.U{7},f_ht.tree.modesizes(v_modes(3)),1),1);
    
    mid = N_mesh(v_modes(1))/2;

    %% Region I: nonnegative/zero modes in vx, vy, and vz
    tree0 = f_ht.tree;
    tree0.modesizes(v_modes(1)) = tree0.modesizes(v_modes(1))/2+1;
    tree0.modesizes(v_modes(2)) = tree0.modesizes(v_modes(2))/2+1;
    tree0.modesizes(v_modes(3)) = tree0.modesizes(v_modes(3))/2+1;
    tree0.freedom(9) = tree0.modesizes(v_modes(1));
    tree0.freedom(7) = tree0.modesizes(v_modes(3));
    tree0.freedom(11) = tree0.modesizes(v_modes(2));
    tree0.freedom(4) = tree0.freedom(8)*tree0.freedom(9);
    tree0.freedom(3) = tree0.freedom(6)*tree0.freedom(7);
    tree0.freedom(5) = tree0.freedom(10)*tree0.freedom(11);
    tree0.freedom(2) = tree0.freedom(4)*tree0.freedom(5);
    tree0.freedom(1) = tree0.freedom(2)*tree0.freedom(3);
    
    % Evaluate the Wigner phase factor on this Fourier block.
    fnc = @(i,j,p,q,l,m) FL_HT3D3V_handle_vec(i,j,p,q,l,m,f_ht,phi_ht,phix_ht,phi_cross,x_extend,k_mod_vx,k_mod_vy,k_mod_vz,h_mesh,N_mesh,n_ghost,dt,H);
    max_rank = tree0.freedom; max_rank(1) = 1;
    % Apply HTACA to construct this Fourier block in HT format.
    f_ht_new = ht.HTACA(fnc,tree0,base_rel_eps,tol_ratio,min_rank,max_rank);
    f_ht_new.U{9}(mid+2:N_mesh(v_modes(1)),:) = 0;
    f_ht_new.U{7}(mid+2:N_mesh(v_modes(3)),:) = 0;
    f_ht_new.U{11}(mid+2:N_mesh(v_modes(2)),:) = 0;
    f_ht_new.tree = f_ht.tree;

    f_ht_new0 = f_ht_new;

    %% Region II: positive vx, negative vy, nonnegative/zero vz
    tree0 = f_ht.tree;
    tree0.modesizes(v_modes(1)) = tree0.modesizes(v_modes(1))/2-1; % Positive vx 
    tree0.modesizes(v_modes(2)) = tree0.modesizes(v_modes(2))/2-1; % Negative vy
    tree0.modesizes(v_modes(3)) = tree0.modesizes(v_modes(3))/2+1; % Neg/Zero vz 
    
    tree0.freedom(9) = tree0.modesizes(v_modes(1));
    tree0.freedom(7) = tree0.modesizes(v_modes(3));
    tree0.freedom(11) = tree0.modesizes(v_modes(2));
    tree0.freedom(4) = tree0.freedom(8)*tree0.freedom(9);
    tree0.freedom(3) = tree0.freedom(6)*tree0.freedom(7);
    tree0.freedom(5) = tree0.freedom(10)*tree0.freedom(11);
    tree0.freedom(2) = tree0.freedom(4)*tree0.freedom(5);
    tree0.freedom(1) = tree0.freedom(2)*tree0.freedom(3);
    
    % Shift j (vx) to positive, shift q (vy) to negative (q+1), leave m (vz) unshifted
    fnc = @(i,j,p,q,l,m) FL_HT3D3V_handle_vec(i,j+mid+1,p,q+1,l,m,f_ht,phi_ht,phix_ht,phi_cross,x_extend,k_mod_vx,k_mod_vy,k_mod_vz,h_mesh,N_mesh,n_ghost,dt,H);
    max_rank = tree0.freedom; max_rank(1) = 1;
    f_ht_new = ht.HTACA(fnc,tree0,base_rel_eps,tol_ratio,min_rank,max_rank);
    
    % vx shift
    f_ht_new.U{9}(mid+2:N_mesh(v_modes(1)),:) = f_ht_new.U{9}(1:mid-1,:);
    f_ht_new.U{9}(1:mid+1,:) = 0;
    
    % vy shift 
    f_ht_new.U{11}(2:mid,:) = f_ht_new.U{11}(1:mid-1,:);
    f_ht_new.U{11}(1,:) = 0;
    f_ht_new.U{11}(mid+1:N_mesh(v_modes(2)),:) = 0;

    % vz padding (was evaluated 1:mid+1)
    f_ht_new.U{7}(mid+2:N_mesh(v_modes(3)),:) = 0;

    f_ht_new.tree = f_ht.tree;
    f_ht_new0 = ht.add(f_ht_new,f_ht_new0); 

    %% Region III: nonnegative/zero vx, positive vy, negative vz
    tree0 = f_ht.tree;
    tree0.modesizes(v_modes(1)) = tree0.modesizes(v_modes(1))/2+1; % Neg/Zero vx
    tree0.modesizes(v_modes(2)) = tree0.modesizes(v_modes(2))/2-1; % Positive vy
    tree0.modesizes(v_modes(3)) = tree0.modesizes(v_modes(3))/2-1; % Negative vz
    
    tree0.freedom(9) = tree0.modesizes(v_modes(1));
    tree0.freedom(7) = tree0.modesizes(v_modes(3));
    tree0.freedom(11) = tree0.modesizes(v_modes(2));
    tree0.freedom(4) = tree0.freedom(8)*tree0.freedom(9);
    tree0.freedom(3) = tree0.freedom(6)*tree0.freedom(7);
    tree0.freedom(5) = tree0.freedom(10)*tree0.freedom(11);
    tree0.freedom(2) = tree0.freedom(4)*tree0.freedom(5);
    tree0.freedom(1) = tree0.freedom(2)*tree0.freedom(3);
    
    % Leave j (vx) unshifted, shift q (vy) to positive, shift m (vz) to negative (m+1)
    fnc = @(i,j,p,q,l,m) FL_HT3D3V_handle_vec(i,j,p,q+mid+1,l,m+1,f_ht,phi_ht,phix_ht,phi_cross,x_extend,k_mod_vx,k_mod_vy,k_mod_vz,h_mesh,N_mesh,n_ghost,dt,H);
    
    max_rank = tree0.freedom; max_rank(1) = 1;
    f_ht_new = ht.HTACA(fnc,tree0,base_rel_eps,tol_ratio,min_rank,max_rank);
    
    % vy shift
    f_ht_new.U{11}(mid+2:N_mesh(v_modes(2)),:) = f_ht_new.U{11}(1:mid-1,:);
    f_ht_new.U{11}(1:mid+1,:) = 0;

    % vz shift
    f_ht_new.U{7}(2:mid,:) = f_ht_new.U{7}(1:mid-1,:);
    f_ht_new.U{7}(1,:) = 0;
    f_ht_new.U{7}(mid+1:N_mesh(v_modes(3)),:) = 0;
    
    % vx padding (was evaluated 1:mid+1)
    f_ht_new.U{9}(mid+2:N_mesh(v_modes(1)),:) = 0;

    f_ht_new.tree = f_ht.tree;
    f_ht_new0 = ht.add(f_ht_new,f_ht_new0);

    %% Region IV: negative vx, nonnegative/zero vy, positive vz
    tree0 = f_ht.tree;
    tree0.modesizes(v_modes(1)) = tree0.modesizes(v_modes(1))/2-1; % Negative vx
    tree0.modesizes(v_modes(2)) = tree0.modesizes(v_modes(2))/2+1; % Neg/Zero vy
    tree0.modesizes(v_modes(3)) = tree0.modesizes(v_modes(3))/2-1; % Positive vz
    
    tree0.freedom(9) = tree0.modesizes(v_modes(1));
    tree0.freedom(7) = tree0.modesizes(v_modes(3));
    tree0.freedom(11) = tree0.modesizes(v_modes(2));
    tree0.freedom(4) = tree0.freedom(8)*tree0.freedom(9);
    tree0.freedom(3) = tree0.freedom(6)*tree0.freedom(7);
    tree0.freedom(5) = tree0.freedom(10)*tree0.freedom(11);
    tree0.freedom(2) = tree0.freedom(4)*tree0.freedom(5);
    tree0.freedom(1) = tree0.freedom(2)*tree0.freedom(3);
    
    % Shift j (vx) to negative (j+1), leave q (vy) unshifted, shift m (vz) to positive
    fnc = @(i,j,p,q,l,m) FL_HT3D3V_handle_vec(i,j+1,p,q,l,m+mid+1,f_ht,phi_ht,phix_ht,phi_cross,x_extend,k_mod_vx,k_mod_vy,k_mod_vz,h_mesh,N_mesh,n_ghost,dt,H);
    max_rank = tree0.freedom; max_rank(1) = 1;
    f_ht_new = ht.HTACA(fnc,tree0,base_rel_eps,tol_ratio,min_rank,max_rank);
    
    % vz shift
    f_ht_new.U{7}(mid+2:N_mesh(v_modes(3)),:) = f_ht_new.U{7}(1:mid-1,:);
    f_ht_new.U{7}(1:mid+1,:) = 0;

    % vx shift
    f_ht_new.U{9}(2:mid,:) = f_ht_new.U{9}(1:mid-1,:);
    f_ht_new.U{9}(1,:) = 0;
    f_ht_new.U{9}(mid+1:N_mesh(v_modes(1)),:) = 0;
    
    % vy padding (was evaluated 1:mid+1)
    f_ht_new.U{11}(mid+2:N_mesh(v_modes(2)),:) = 0;

    f_ht_new.tree = f_ht.tree;
    f_ht_new0 = ht.add(f_ht_new,f_ht_new0);
    
    % Use conjugate symmetry to fill the complementary Fourier modes.
    f_ht_new1 = f_ht_new0;
    f_ht_new1.B{1} = conj(f_ht_new1.B{1});
    f_ht_new1.B{2} = conj(f_ht_new1.B{2});
    f_ht_new1.B{3} = conj(f_ht_new1.B{3});
    f_ht_new1.B{4} = conj(f_ht_new1.B{4});
    f_ht_new1.B{5} = conj(f_ht_new1.B{5});

    f_ht_new1.U{8} = conj(f_ht_new1.U{8});
    f_ht_new1.U{6} = conj(f_ht_new1.U{6});
    f_ht_new1.U{10} = conj(f_ht_new1.U{10});

    f_ht_new1.U{9} = conj(f_ht_new1.U{9}(opp_index_vx(1:N_mesh(v_modes(1))),:));
    f_ht_new1.U{7} = conj(f_ht_new1.U{7}(opp_index_vz(1:N_mesh(v_modes(3))),:));
    f_ht_new1.U{11} = conj(f_ht_new1.U{11}(opp_index_vy(1:N_mesh(v_modes(2))),:));
 
    f_ht_new0 = ht.add(f_ht_new0,f_ht_new1);

    %% Remove the duplicated zero-frequency contribution
    f_ht_new = f_ht_new0;
    f_ht_new.U{9}(2:mid,:) = 0;
    f_ht_new.U{9}(mid+2:end,:) = 0;
    f_ht_new.U{7}(2:mid,:) = 0;
    f_ht_new.U{7}(mid+2:end,:) = 0;
    f_ht_new.U{11}(2:mid,:) = 0;
    f_ht_new.U{11}(mid+2:end,:) = 0;
    f_ht_new = ht.multiply_scalar(f_ht_new,-1);

    f_ht_new0 = ht.add(f_ht_new0,f_ht_new);

    %% Add back the zero-frequency velocity mode
    f_ht_new = f_ht;
    f_ht_new.U{9}(1:mid,:) = 0;
    f_ht_new.U{9}(mid+2:end,:) = 0;
    f_ht_new.U{7}(1:mid,:) = 0;
    f_ht_new.U{7}(mid+2:end,:) = 0;
    f_ht_new.U{11}(1:mid,:) = 0;
    f_ht_new.U{11}(mid+2:end,:) = 0;

    f_ht_new = ht.add(f_ht_new0,f_ht_new);

    f_ht_new.U{11} = ifft(ifftshift(f_ht_new.U{11},1),f_ht_new.tree.modesizes(v_modes(2)),1);
    f_ht_new.U{7} = ifft(ifftshift(f_ht_new.U{7},1),f_ht_new.tree.modesizes(v_modes(3)),1);
    f_ht_new.U{9} = ifft(ifftshift(f_ht_new.U{9},1),f_ht_new.tree.modesizes(v_modes(1)),1);
end
