function rho = compute_rho(f_ht,v_modes, h_mesh)
%COMPUTE_RHO Compute the 2D charge density from a 2D2V HT distribution.
%
%   rho = compute_rho(f_ht, v_modes, h_mesh)
%
%   This function integrates out the velocity variables from f(x,v) using
%   the midpoint rule. The remaining modes are the spatial dimensions.
%
%   Inputs:
%       f_ht     - 4D HTD structure representing the phase-space distribution
%       v_modes  - row vector specifying the indices of velocity dimensions 
%                  in f_ht.tree
%       h_mesh   - grid spacing in each velocity dimension
%
%   Output:
%       rho - 2D HTD structure after contracting over velocity dimensions


    % Set function handles for contracting velocity modes using the midpoint rule.
    D = length(v_modes);    % number of velocity dimensions
    F_list = cell(1,D);
    for d = 1: D
        F_list{d} = @(x) h_mesh(d)*sum(x,1);  % midpoint rule integration
        % F_list{d} = @(x) h_mesh(d)*(sum(x,1)-0.5*(x(1,:)-x(end,:))); 
    end
    
    % Apply function handles to specified modes and squeeze out reduced dimensions.
    rho = ht.squeeze(ht.mode_apply_function(f_ht,F_list,v_modes));

end
