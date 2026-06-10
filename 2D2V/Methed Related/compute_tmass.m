function total_mass = compute_tmass(f_ht,h_mesh)
%COMPUTE_TMASS Compute the total mass (integral) of an HT-represented scalar field
%
%   total_mass = compute_tmass(f_ht, h_mesh)
%
%   This function computes the total mass of a scalar function f(x₁, ..., x_D)
%   represented in Hierarchical Tucker (HT) format over a uniform grid.
%   The integral is approximated using the midpoint rule in each dimension.
%
%   Specifically, the total mass is computed as:
%       ∫ f(x₁, ..., x_D) dx ≈ ∑_i f_i × (∏_d h_mesh(d))
%   where h_mesh is the grid spacing in each dimension.
%
%   INPUTS:
%     f_ht    - HTD structure of the scalar function f with fields .U, .B, .tree, .rank
%     h_mesh  - grid spacing in each spatial dimension (1×D vector)
%
%   OUTPUT:
%     total_mass - scalar value representing the total mass ∫ f(x) dx

    % Number of physical dimensions
    D = f_ht.tree.orders(1);
    
    % Define integration function handles for each mode (midpoint rule)
    F_list = cell(1, D);
    index_ = cell(1, D);
    for d = 1:D
        F_list{d}  = @(x) h_mesh(d) * sum(x, 1);  % apply 1D midpoint integration
        index_{d}  = 1;                           % scalar index to evaluate final result
    end

    % Apply integration along all modes, then evaluate scalar at index (1,1,...)
    total_mass = ht.evaluate_index(ht.mode_apply_function(f_ht,F_list,1:D),index_{:});

end