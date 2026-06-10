function [dom_min, dom_max, dom_len, h_mesh, x, v, x_extend, v_extend, kx] = set_spatial_grid_info_6d(N_mesh, n_ghost, x_modes, v_modes, prob)
%SET_SPATIAL_GRID_INFO_6D Set grid information for the 3D3V phase-space domain
%
%   [dom_min, dom_max, dom_len, h_mesh, x, v, x_extend, v_extend, kx]
%       = set_spatial_grid_info_6d(N_mesh, n_ghost, x_modes, v_modes, prob)
%
%   This function constructs spatial and velocity grids, ghost-cell extensions,
%   and FFT-compatible frequency grids for simulations in 6D phase space.
%
%   Inputs:
%       N_mesh  - A vector specifying the number of grid points in each dimension
%       n_ghost - Number of ghost cells used for boundary extensions
%       x_modes - Vector of mode indices corresponding to spatial dimensions
%       v_modes - Vector of mode indices corresponding to velocity dimensions 
%       prob    - Problem type identifier
%                  1 or 2: k = 0.5, v_max = 2*pi
%                  3:      k = 0.2, v_max = 8
%                  4:      k = 0.5, v_max = 2*pi
%                  5:      k = 0.3, v_max = 8
%
%   Outputs:
%       dom_min  - Minimum coordinates of the computational domain
%       dom_max  - Maximum coordinates of the computational domain
%       dom_len  - Domain size in each dimension
%       h_mesh   - Grid spacing in each dimension
%       x        - 1 x 3 cell arrays containing 
%                  1D arrays of spatial grid points in x, y and z directions
%       v        - 1 x 3 cell arrays containing 
%                  1D arrays of velocity grid points in vx, vy and vz directions
%       x_extend - 1 x 3 cell arrays containing
%                  extended grid with ghost cells in x, y and z directions
%       v_extend - 1 x 3 cell arrays containing
%                  extended grid with ghost cells in vx, vy and vz directions
%       kx       - 1 x 3 cell arrays containing
%                    fundamental frequencies in x, y and z directions
%
%   Note:
%       All grid points are centered at cell centers.
%       Extended arrays contain [grid_coordinate, logical_index], where ghost zones
%       are assigned boundary indices or wrapped indices depending on the type (spatial/velocity)

% Set problem-specific parameters
if prob == 1 || prob == 2
    k    = 0.5;      % Wave number or scaling parameter for spatial domain
    vmax = 2*pi;     % Maximum velocity in velocity directions
elseif prob == 3
    k    = 0.2;
    vmax = 8;
elseif prob == 4
    k    = 0.5;
    vmax = 2*pi;
elseif prob == 5
    k    = 0.3;
    vmax = 8;
end

% Set the total number of dimensions
D = 6;

% Define computational domain bounds
dom_min = zeros(1,D);
dom_max = zeros(1,D);
dom_min(x_modes) = 0;
dom_max(x_modes) = 2*pi/k;
dom_min(v_modes) = -vmax;
dom_max(v_modes) = vmax;

% Compute domain length and grid spacing
dom_len = dom_max - dom_min;  % Total length in each direction
h_mesh  = dom_len ./ N_mesh;  % Grid step size in each dimension

% Generate uniform grid points (centered at cell centers)
x = cell(1,3); v = cell(1,3);
for d = 1: 3
    x{d} = h_mesh(x_modes(d)) * (0:N_mesh(x_modes(d))-1)' + dom_min(x_modes(d));
    v{d} = h_mesh(v_modes(d)) * (0:N_mesh(v_modes(d))-1)' + dom_min(v_modes(d));
end

% Initialize extended grid arrays for ghost cells
x_extend = cell(1,3); v_extend = cell(1,3);
for d = 1: 3
    x_extend{d} = zeros(N_mesh(x_modes(d)) + 2 * n_ghost, 2);
    v_extend{d} = zeros(N_mesh(v_modes(d)) + 2 * n_ghost, 2);
end

% Extend spatial and velocity grids with ghost cells
for d = 1: 3
    x_extend{d}(:,1) = h_mesh(x_modes(d)) * (-n_ghost : N_mesh(x_modes(d))-1+n_ghost)'  + dom_min(x_modes(d));    
    v_extend{d}(:,1) = h_mesh(v_modes(d)) * (-n_ghost : N_mesh(v_modes(d))-1+n_ghost)' + dom_min(v_modes(d));   
end

% Fill in interior indices (physical domain)
for d = 1: 3
    x_extend{d}(n_ghost+1 : N_mesh(x_modes(d))+n_ghost, 2)  = (1:N_mesh(x_modes(d)))';
    v_extend{d}(n_ghost+1 : N_mesh(v_modes(d))+n_ghost, 2)  = (1:N_mesh(v_modes(d)))';
end

% Fill in right ghost cells (beyond domain boundary)
for d = 1: 3
    x_extend{d}(N_mesh(x_modes(d)) + n_ghost + 1 : end, 2) = (1:n_ghost)';
    v_extend{d}(N_mesh(v_modes(d)) + n_ghost + 1 : end, 2) = N_mesh(v_modes(d));  % Last velocity index repeated
end

% Fill in left ghost cells (before domain boundary)
for d = 1: 3
    x_extend{d}(1:n_ghost, 2) = (N_mesh(x_modes(d)) - n_ghost + 1 : N_mesh(x_modes(d)))';
    v_extend{d}(1:n_ghost, 2) = 1;  % First velocity index repeated
end

% Create the arrays of fundamental frequencies
kx = cell(1,3);
for d = 1: 3
    kx{d} = (2*pi/dom_len(x_modes(d)))*fftfreq(N_mesh(x_modes(d))); 
end

end

function [k] = fftfreq(n)
% Returns the Discrete Fourier Transform sample frequencies.
% The returned frequencies can be scaled and used for differentiation.
%
% We only need to consider two cases: (1) n is even and (2) n is odd.

    if mod(n,2) == 0
        k = [0:(n/2-1) -n/2:-1]';
    else
        k = [0:(n-1)/2 -(n-1)/2:-1]';
    end

end
