% EXAMPLE_TRUNCATE_AND_GRAMIAN.M
% Demonstrates how to orthogonalize, compute Gramians, and truncate
% a Hierarchical Tucker decomposition using ht.orthogonalize, 
% ht.gramians_orthog, and ht.truncate. This example uses a compressible 
% 3D3V Maxwellian function as a benchmark.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

%% === 1. Build the dimension tree and define a compressible HTD ===

% Define moderate mode sizes for each of the 6 dimensions (3 spatial + 3 velocity)
modesizes = [32, 64, 48, 32, 64, 16];

% Construct a balanced binary tree with 6 leaves and 11 total nodes
children = [2,3; 4,5; 6,7; 8,9; 10,11; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0];
dim2ind = [8,9,10,11,6,7];  % mapping of physical dimensions to leaf nodes
tree = ht.make_tree(children, dim2ind, modesizes);

% Define structured 3D spatial and velocity grids
hx  = 4*pi/modesizes(1); hy  = 4*pi/modesizes(2);  hz  = 4*pi/modesizes(3);
hv1 = 4*pi/modesizes(4); hv2 = 4*pi/modesizes(5);  hv3 = 4*pi/modesizes(6);
x  = ((1:modesizes(1))*hx - 0.5*hx)';
y  = ((1:modesizes(2))*hy - 0.5*hy)';
z  = ((1:modesizes(3))*hz - 0.5*hz)';
v1 = ((1:modesizes(4))*hv1 - 0.5*hv1 - 2*pi)';
v2 = ((1:modesizes(5))*hv2 - 0.5*hv2 - 2*pi)';
v3 = ((1:modesizes(6))*hv3 - 0.5*hv3 - 2*pi)';

% Define 3D3V BGK Maxwellian as a function handle
maxwellian_bgk = @(i,j,p,q,m,n) local_maxwellian(x(i),y(j),z(p),v1(q),v2(m),v3(n));

% Run HTACA to compress the Maxwellian function in HT format
base_rel_tol = 1e-6;
tol_ratio = 0.1;
min_rank = ht.get_default_min_rank(tree);
max_rank = ht.get_default_max_rank(tree);
A = ht.HTACA(maxwellian_bgk, tree, base_rel_tol, tol_ratio, min_rank, max_rank);

disp('--- Original HTD ---');
disp(['Rank: ', mat2str(A.rank)]);
disp(['Frobenius norm: ', num2str(ht.norm(A))]);

%% === 2. Orthogonalize the HTD ===

A_orth = ht.orthogonalize(A);

%% === 3. Compute Gramians for orthogonalized HTD ===

G = ht.gramians_orthog(A_orth);
disp(['Computed Gramians for ', num2str(length(G)), ' nodes.']);

%% === 4. Truncate HTD with a prescribed tolerance ===

tol = 1e-2;  % truncation tolerance (relative to original norm)
A_trunc = ht.truncate(A, tol * ht.norm(A), min_rank, max_rank);

% Compute relative error of truncation
A_error = ht.add(A, ht.multiply_scalar(A_trunc, -1));

disp('--- Truncated HTD ---');
disp(['New rank: ', mat2str(A_trunc.rank)]);
disp(['Relative error: ', num2str(ht.norm(A_error) / ht.norm(A))]);


%% === Local Maxwellian Function ===
function val = local_maxwellian(x, y, z, vx, vy, vz)
    % Local Maxwellian with space-dependent density and velocity, constant temperature
    rho = 1 + 0.2 * sin(2*pi*x) .* sin(2*pi*y) .* sin(2*pi*z);  % density
    u1  = 0.1 * sin(pi*x);  % velocity field in x
    u2  = 0.1 * sin(pi*y);  % velocity field in y
    u3  = 0.1 * sin(pi*z);  % velocity field in z
    T   = 1;                % constant temperature

    % Compute squared relative velocity
    v_rel_sq = (vx - u1).^2 + (vy - u2).^2 + (vz - u3).^2;

    % Evaluate Maxwellian
    val = rho ./ ( (2*pi*T)^(3/2) ) .* exp( -v_rel_sq ./ (2*T) );
end
