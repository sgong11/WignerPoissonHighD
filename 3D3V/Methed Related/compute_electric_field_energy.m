function E_energy = compute_electric_field_energy(h_mesh, E_field)
%COMPUTE_ELECTRIC_FIELD_ENERGY Compute the total electric field energy from 
% HT-represented electric field components.
%
%   E_energy = compute_electric_field_energy(h_mesh, E_field)
%
%   This function calculates the electric field energy using the formula:
%       E_energy = ½ ∫ (|E₁|² + |E₂|² (+ |E₃|²)) dx
%   where the electric field components are represented as HTD tensors.
%
%   Inputs:
%     h_mesh  - grid spacing in each spatial dimension (1×D vector)
%     E_field - 1 x 2 or 1 x 3 cell array containing electric field
%               components in HT format
%   Output:
%     E_energy - scalar value of total electric field energy

    % Determine the number of spatial dimensions: D = 2 or 3
    D = length(h_mesh);

    % Initialize energy accumulator
    E_energy = 0;

    % Accumulate energy from each component: ||E_d||²
    for d = 1: D
        % ht.norm returns the Frobenius norm of the HT tensor
        E_energy = E_energy + ht.norm(E_field{d})^2;
    end

    % Multiply by 1/2 * dx * dy * (dz) to get total energy.
    E_energy = 0.5*prod(h_mesh)*E_energy;

end
