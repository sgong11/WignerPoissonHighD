function energy_kinetic = compute_kinetic_energy(f_ht,h_mesh,x_modes,v_modes,v)
%COMPUTE_KINETIC_ENERGY Compute the kinetic energy of a high-dimensional distribution f in HT format
%
%   energy_kinetic = compute_kinetic_energy(f_ht,h_mesh,x_modes,v_modes,v)
%
%   This function computes the total kinetic energy by integrating:
%
%       energy_kinetic = ½ ∫∫ f(x,v) · (vₓ² + vᵧ² (+ v_z²)) dx dv
%
%   where f is a phase-space distribution represented in Hierarchical Tucker (HT) format.
%
%   Inputs:
%       f_ht    - 6D HTD structure with fields .U, .B, .tree, .rank
%              representing the phase-space distribution function f
%       h_mesh  - grid spacing in each dimension (1×D vector)
%       x_modes - vector of mode indices corresponding to spatial dimensions
%       v_modes - vector of mode indices corresponding to velocity dimensions
%       v       - 1 x 2 or 1 x 3 cell array containing 1D velocity grid arrays in
%                     vx, vy (and vz) directions
%
%   Outputs:
%        energy - scalar value representing total kinetic energy
    
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
    energy_kinetic_components = zeros(1,length(x_modes));

    for i = 1: length(v_modes)
        d      = v_modes(i);                 % velocity mode index

        % Clone default integration handles and insert weighted velocity integral.
        F_list_d = F_list;
        F_list_d{d} = @(x) h_mesh(d) * sum(x .* v{i}(:).^2, 1);  % apply v_i^2 * f

        % Evaluate integrated scalar: integral f(x,v) v_i^2 dx dv.
        energy_kinetic_components(d) = ht.evaluate_index(ht.mode_apply_function(f_ht,F_list_d,1:D),index_{:});
    end

    % Assign output
    energy_kinetic = 0.5 * sum(energy_kinetic_components);

end
