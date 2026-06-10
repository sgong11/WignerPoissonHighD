function momentum = compute_momentum(f_ht,h_mesh,x_modes,v_modes,v)
%COMPUTE_MOMENTUM Compute the momentum components from a high-dimensional distribution f in HT format
%
%   momentum = compute_momentum(f_ht,h_mesh,x_modes,v_modes,v)
%
%   This function computes the momentum components by integrating the
%   product f(x,v) · v_i over all spatial and velocity variables, i.e.:
%
%       momentum_i = ∫∫ f(x,v) · v_i(x) dx dv
%
%   Inputs:
%       f_ht    - 6D HTD structure with fields .U, .B, .tree, .rank
%              representing the phase-space distribution function f
%       h_mesh  - grid spacing in each dimension (1×D vector)
%       x_modes - vector of mode indices corresponding to spatial dimensions
%       v_modes - vector of mode indices corresponding to velocity dimensions
%       v       - 1 x 2 or 1 x 3 cell array containing 1D velocity grid arrays 
%                   in vx, vy (and vz) directions
%
%   Outputs:
%       momentum - 1 x 2 or 1 x 3 array containing the 2 or 3 momentum components
    
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

        % Clone default integration handles and insert weighted velocity integral.
        F_list_d = F_list;
        F_list_d{d} = @(x) h_mesh(d) * sum(x .* v{i}(:), 1);  % apply v_i * f

        % Evaluate integrated scalar: integral f(x,v) v_i dx dv.
        momentum(i) = ht.evaluate_index(ht.mode_apply_function(f_ht,F_list_d,1:D),index_{:});
    end

end
