function E_field = get_E_field_FFT(rho,dom_len,h_mesh,base_rel_tol, tol_ratio, min_rank, max_rank, kx)
%GET_E_FIELD_FFT Compute the periodic electric field in HT format.
%
%   E_field = get_E_field_FFT(rho,dom_len,h_mesh,base_rel_tol, tol_ratio, min_rank, max_rank, kx)
% 
% Inputs:
%   rho         - 2D HTD charge density structure
%   dom_len     - domain size in each dimension (length D vector)
%   h_mesh      - grid spacing in each dimension (length D vector)
%   base_rel_tol - baseline relative truncation tolerance for HTACA
%   tol_ratio   - relative decay factor for HTACA tolerance across levels
%   min_rank    - minimum allowed HT rank for each node (vector of length N_node)
%   max_rank    - maximum allowed HT rank for each node (vector of length N_node)
%   kx          - 1 x 2 cell array of Fourier frequency vectors in x and y
%
% Output:
%   E_field - 1 x 2 cell array of electric field components in HT format

% Note: For periodic problems, we have a solvability condition 
% requiring the domain average of rho to be zero. This is enforced by subtracting 
% the mean value μ:
%     rho := rho - μ,
% where μ is the domain average computed using total mass and grid spacing.

    % Enforce solvability condition.
    mu  = 1./prod(dom_len)*compute_tmass(rho,h_mesh);
    rho = ht.add_scalar(rho,-mu);
    
    % Apply FFT to the bases of \rho
    tree = rho.tree;
    D    = tree.orders(1);
    for d = 1: D
        node        = tree.dim2ind(d);
        rho.U{node} = fft(rho.U{node}, tree.modesizes(d), 1);
    end
    
    % Construct the function handle for solving Δφ = rho in the frequency domain.
    if D == 2
        data = @(i,j) phi_hat_fnc(rho,kx,D,i,j); % 2D case
    else
        data = @(i,j,k) phi_hat_fnc(rho,kx,D,i,j,k); % Higher-dimensional fallback
    end
    
    % Construct φ̂  in HT format using the function handle and prescribed
    % truncation settings (via HTACA)
    phi_hat = ht.HTACA(data,tree,base_rel_tol, tol_ratio, min_rank, max_rank);

    % Apply inverse FFT to the bases of \hat{\phi} and recover \phi
    phi = phi_hat;
    for d = 1: D
        node = tree.dim2ind(d);
        phi.U{node} = ifft(phi_hat.U{node}, tree.modesizes(d), 1);
    end

    % Compute E = -grad(phi) as the inverse FFT of -i*k*phi_hat.
    E_field = cell(1,D);
    E_field = cellfun(@(x) phi, E_field, 'UniformOutput', false);
    for d = 1: D
        node = tree.dim2ind(d);
        E_field{d}.U{node} = ifft(phi_hat.U{node}.*(-1i*kx{d}), tree.modesizes(d), 1);
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
