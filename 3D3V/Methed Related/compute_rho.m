function rho = compute_rho(f_ht,v_modes, h_mesh)
%COMPUTE_RHO Compute the charge density ρ from a high-dimensional distribution f in HT format
%
%   rho = compute_rho(f_ht, v_modes, h_mesh)
%
%   This function computes the charge density ρ by integrating out the velocity
%   variables from the distribution function f using the midpoint rule.
%   The input f_ht is assumed to be a high-dimensional tensor in HT format,
%   where some dimensions correspond to physical space (e.g., x, y, z) and others
%   to velocity space (e.g., vx, vy, vz).
%
%   Inputs:
%       f_ht     - 6D HTD structure with fields .U, .B, .tree, .rank
%                  representing the phase-space distribution function f
%       v_modes  - row vector specifying the indices of velocity dimensions 
%                  in f_ht.tree
%       h_mesh   - grid spacing in each velocity dimension (1×D_v vector)
%
%   Output:
%       rho - 2D or 3D HTD structure after contracting over velocity dimensions;
%             retains only spatial dimensions (x, y, z in 3D3V)


    % Set function handles for contracting velocity modes with the midpoint rule.
    D = length(v_modes);    % number of velocity dimensions
    F_list = cell(1,D);
    for d = 1: D
        F_list{d} = @(x) h_mesh(d)*sum(x,1);  % midpoint rule integration
    end
    
    % Apply function handles to specified modes and squeeze out reduced dimensions.
    rho = ht.squeeze(ht.mode_apply_function(f_ht,F_list,v_modes));

end
