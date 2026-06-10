function [phi,phix_ht,phixy_ht] = get_phi_FFT(rho,base_rel_tol, tol_ratio, min_rank, max_rank, kx)
%GET_PHI_FFT Solve the periodic Poisson equation in HT format.
%
%   [phi, phix_ht, phixy_ht] =
%       get_phi_FFT(rho, base_rel_tol, tol_ratio, min_rank, max_rank, kx)
% 
% Inputs:
%   rho         - 2D HTD charge density structure
%   base_rel_tol - baseline relative truncation tolerance for HTACA
%   tol_ratio   - relative decay factor for HTACA tolerance across levels
%   min_rank    - minimum allowed HT rank for each node (vector of length N_node)
%   max_rank    - maximum allowed HT rank for each node (vector of length N_node)
%   kx          - 1 x 2 cell array of Fourier frequency vectors in x and y
%
% Outputs:
%   phi      - electrostatic potential in HT format
%   phix_ht  - 1 x 2 cell array of first derivatives of phi
%   phixy_ht - mixed derivative of phi used by clamped interpolation

% The drivers normalize mass to the spatial domain volume, so the neutral
% background has density one. Subtracting one enforces the periodic Poisson
% solvability condition.

    % Enforce solvability condition.
    rho = ht.add_scalar(rho,-1);
    
    % Apply FFT to the bases of \rho
    tree = rho.tree;
    D    = tree.orders(1);
    for d = 1: D
        node        = tree.dim2ind(d);
        rho.U{node} = fft(rho.U{node}, tree.modesizes(d), 1);
    end
    
    % Construct the function handle for solving Δφ = ρ in the frequency domain.
    if D == 2
        data = @(i,j) phi_hat_fnc(rho,kx,D,i,j); % 2D case
    else
        data = @(i,j,k) phi_hat_fnc(rho,kx,D,i,j,k); % Higher-dimensional fallback
    end
    

    % Construct φ̂  in HT format using the function handle and prescribed
    % truncation settings (via HTACA)
    phi_hat = ht.HTACA(data,tree,base_rel_tol, tol_ratio, min_rank, max_rank);
    % 
    % % Apply inverse FFT to the bases of \hat{\phi} and recover \phi
    phi = phi_hat;
    for d = 1: D
        node = tree.dim2ind(d);
        phi.U{node} = ifft(phi_hat.U{node}, tree.modesizes(d), 1);
    end

    % Compute first derivatives as the inverse FFT of i*k*phi_hat.
    phix_ht = cell(1,D);
    phix_ht = cellfun(@(x) phi, phix_ht, 'UniformOutput', false);
    for d = 1: D
        node = tree.dim2ind(d);
        phix_ht{d}.U{node} = ifft(phi_hat.U{node}.*(1i*kx{d}), tree.modesizes(d), 1);
    end
    
    phixy_ht = phi_hat;
    for d = 1: D
        node = tree.dim2ind(d);
        phixy_ht.U{node} = ifft(phi_hat.U{node}.*(1i*kx{d}), tree.modesizes(d), 1);
    end
    

end

function value = phi_hat_fnc(rho,k_fre,D,varargin)
%PHI_HAT_FNC Compute φ̂ = ρ̂ / (k₁² + k₂² + ⋯) in frequency space for Poisson solver.
%
%   This function evaluates the right-hand side in Fourier space:
%     φ̂(k) = ρ̂(k) / (k₁² + k₂² + ⋯)
%   with zero padding at k = 0 to enforce solvability.

    index_len = [];
    K = 0;
    for d = 1: D
        % Verify the sizes of the input indices
        index_len = [index_len length(varargin{d})];

        % Form the denominator of 1./(kx^2 + ky^2 + kz^2)
        K = K + (k_fre{d}(varargin{d})).^2;
    end
    N = max(index_len);

    % Form 1./(kx^2 + ky^2 + kz^2) (when kx=ky=kz=0, this is set as 0)
    one_over_K = zeros(N,1);
    one_over_K(K>0) = 1./K(K>0);

    if sum(index_len > 1) == 1
        % Case 1: input indices are long a fiber
        value = one_over_K.*ht.evaluate_fiber(rho,varargin{:});
        return
    else
        % Case 2: input indices are not along a fiber
        varargin = cellfun(@(x) x.*ones(N,1), varargin, 'UniformOutput',false);
        value = one_over_K.*ht.evaluate_index(rho,varargin{:});
    end
end
