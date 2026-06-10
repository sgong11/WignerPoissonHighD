function E_field = get_E_field_FFT(rho,dom_len,h_mesh,base_rel_tol, tol_ratio, min_rank, max_rank, kx)
%GET_E_FIELD_FFT Compute the electric field E1, E2 (and E3 for 3D) by solving the Poisson
% equation -Δφ = ρ and computing E = -∇φ using the FFT.
%
%   E_field = get_E_field_FFT(rho,dom_len,h_mesh,base_rel_tol, tol_ratio, min_rank, max_rank, kx)
% 
% Inputs:
%   rho         - 2D or 3D HTD structure with fields .U, .B, .tree, .rank
%   dom_len     - domain size in each dimension (length D vector)
%   h_mesh      - grid spacing in each dimension (length D vector)
%   base_rel_tol - baseline relative truncation tolerance for HTACA
%   tol_ratio   - relative decay factor for HTACA tolerance across levels
%   min_rank    - minimum allowed HT rank for each node (vector of length N_node)
%   max_rank    - maximum allowed HT rank for each node (vector of length N_node)
%   kx          - 1 x 2 or 1 x 3 cell arrays containing
%                    fundamental frequencies in x, y (and z) directions
%
% Output:
%   E_field - 1 x 2 or 1 x 3 cell array containing electric field components
%             in HT format
%
% Note: For periodic problems, we have a solvability condition 
% requiring the domain average of ρ to be zero. This is enforced by subtracting 
% the mean value μ:
%     ρ := ρ - μ,
% where μ is the domain average of ρ computed using total mass and grid spacing.

    % Enforce the solvability condition.
    mu  = 1./prod(dom_len)*compute_tmass(rho,h_mesh);
    rho = ht.add_scalar(rho,-mu);
    
    % Apply FFT to the bases of \rho
    tree = rho.tree;
    D    = tree.orders(1);
    for d = 1: D
        node        = tree.dim2ind(d);
        rho.U{node} = fft(rho.U{node}, tree.modesizes(d), 1);
    end
    
    % Construct the function handle for solving Δφ = ρ in the frequency domain
    if d == 2
        data = @(i,j) phi_hat_fnc(rho,kx,D,i,j); % 2d case
    else
        data = @(i,j,k) phi_hat_fnc(rho,kx,D,i,j,k); % 3d case
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

    % Compute E as the inverse FFT of -i*k*phi_hat in each direction.
    E_field = cell(1,D);
    E_field = cellfun(@(x) phi, E_field, 'UniformOutput', false);
    for d = 1: D
        node = tree.dim2ind(d);
        E_field{d}.U{node} = ifft(phi_hat.U{node}.*(-1i*kx{d}), tree.modesizes(d), 1);
    end

end

function value = phi_hat_fnc(rho,k_fre,D,varargin)
%PHI_HAT_FNC Compute phi_hat = rho_hat / (k1^2 + k2^2 + ...) in frequency space.
%
%   This function evaluates the right-hand side in Fourier space:
%     phi_hat(k) = rho_hat(k) / (k1^2 + k2^2 + ...)
%   with zero padding at k = 0 to enforce solvability.

    index_len = [];
    K = 0;
    for d = 1: D
        % Record the sizes of the input indices.
        index_len = [index_len length(varargin{d})];

        % Form the denominator of 1./(kx^2 + ky^2 + kz^2)
        K = K + (k_fre{d}(varargin{d})).^2;
    end
    N = max(index_len);

    % Form 1./(kx^2 + ky^2 + kz^2) (when kx=ky=kz=0, this is set as 0)
    one_over_K = zeros(N,1);
    one_over_K(K>0) = 1./K(K>0);

    if sum(index_len > 1) == 1
            % Case 1: input indices lie along a fiber.
        value = one_over_K.*ht.evaluate_fiber(rho,varargin{:});
        return
    else
        % Case 2: input indices do not lie along a fiber.
        varargin = cellfun(@(x) x.*ones(N,1), varargin, 'UniformOutput',false);
        value = one_over_K.*ht.evaluate_index(rho,varargin{:});
    end
end
