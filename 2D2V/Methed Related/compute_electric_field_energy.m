function E_energy = compute_electric_field_energy(h_mesh, E_field)
%COMPUTE_ELECTRIC_FIELD_ENERGY Compute the total electric field energy from 
% HT-represented electric field components.
%
%   E_energy = compute_electric_field_energy(h_mesh, E_field)
%
%   This function calculates the electric field energy using the formula:
%       E_energy = 1/2 * integral(|E_1|^2 + |E_2|^2) dx dy
%   where the electric field components are represented as HTD tensors.
%
%   Inputs:
%     h_mesh  - grid spacing in each spatial dimension
%     E_field - 1 x 2 cell array of electric field components in HT format
%   Output:
%     E_energy - scalar value of total electric field energy

    % Determine the number of spatial dimensions.
    D = length(h_mesh);

    % Initialize energy accumulator
    E_energy = 0;

    % Accumulate energy from each component.
    for d = 1: D
        % ht.norm returns the Frobenius norm of the HT tensor
        E_energy = E_energy + ht.norm(E_field{d})^2;
    end

    % Multiply by the spatial volume element and the 1/2 prefactor.
    E_energy = 0.5*prod(h_mesh)*E_energy;

end
