function momentum = compute_momentum(f_ht,h_mesh,x_modes,v_modes,v)
%COMPUTE_MOMENTUM Compute momentum components of a 2D2V HT distribution.
%
%   momentum = compute_momentum(f_ht,h_mesh,x_modes,v_modes,v)
%
%   This function computes the momentum components by integrating the
%   product f(x,v) · v_i over all spatial and velocity variables, i.e.:
%
%       momentum_i = integral f(x,v) * v_i dx dv
%
%   Inputs:
%       f_ht    - 4D HTD structure representing the phase-space distribution
%       h_mesh  - grid spacing in each dimension
%       x_modes - vector of mode indices corresponding to spatial dimensions
%       v_modes - vector of mode indices corresponding to velocity dimensions
%       v       - 1 x 2 cell array containing velocity grids in vx and vy
%
%   Outputs:
%       momentum - 1 x 2 array containing the two momentum components
    
% Total number of physical modes (space + velocity)
    D = f_ht.tree.orders(1);

    % Midpoint integration for each dimension
    F_list = cell(1, D);
    index_ = cell(1, D);
    for d = 1:D
        F_list{d}  = @(x) h_mesh(d) * sum(x, 1);  % apply 1D midpoint integration
        index_{d}  = 1;                           % scalar index to evaluate final result
    end

    % Initialize output container
    momentum = zeros(1,length(x_modes));

    for i = 1: length(v_modes)
        d      = v_modes(i);                 % velocity mode index

        % Clone default integration handles and insert weighted velocity integral
        F_list_d = F_list;
        F_list_d{d} = @(x) h_mesh(d) * sum(x .* v{i}(:), 1);  % vectorized: apply v_i(x)·f

        % Evaluate integrated scalar: ∫∫ f(x,v) v_i dx dv
        momentum(i) = ht.evaluate_index(ht.mode_apply_function(f_ht,F_list_d,1:D),index_{:});
    end

end
