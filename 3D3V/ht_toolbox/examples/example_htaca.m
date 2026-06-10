% EXAMPLE_HTACA.M
% Demonstrates HTACA on a structured 6D Maxwellian and evaluates accuracy
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..')));

%--------------------------------------------------------------------------
% Input Function Handle Specification for HTACA
%
%   The input 'data' must be a function handle representing a d-dimensional
%   tensor with the following signature:
%
%       data = @(i1, i2, ..., id) ...
%
%   where:
%       - i1, i2, ..., id are indices corresponding to each of the d modes.
%       - The function must support two access modes:
%
%         (1) Random access mode:
%             All inputs i1, ..., id are column vectors of the same length N.
%             The function returns a column vector of length N:
%                 output = [data(i1(k), ..., id(k)) for k = 1:N]
%
%         (2) Fiber access mode:
%             One of the indices (say, iμ) is a column vector of length N,
%             while all other indices are scalars. The function should return
%             an N×1 vector:
%                 output = [data(iμ(k), ..., fixed indices) for k = 1:N]
%
%   Output:
%       - A column vector of evaluated values corresponding to the given
%         multi-index query.
%
%   Note:
%       - It is strongly recommended that users optimize the fiber access
%         mode if possible. For many physical models or data sources, fiber
%         evaluations can be made significantly faster than general random
%         access. Efficient fiber support will greatly accelerate HTACA
%         construction.
%--------------------------------------------------------------------------

% === 1. Build dimension tree ===
modesizes = [32, 64, 48, 32, 64, 16];  % 6D: (x,y,z,v1,v2,v3)
children = [2,3; 4,5; 6,7; 8,9; 10,11; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0];
dim2ind = [8,9,10,11,6,7];
tree = ht.make_tree(children, dim2ind, modesizes);

% === 2. Build spatial and velocity grids ===
hx  = 4*pi/modesizes(1); hy  = 4*pi/modesizes(2);  hz  = 4*pi/modesizes(3);
hv1 = 4*pi/modesizes(4); hv2 = 4*pi/modesizes(5);  hv3 = 4*pi/modesizes(6);
x  = ((1:modesizes(1))*hx - 0.5*hx)';     % x-grid
y  = ((1:modesizes(2))*hy - 0.5*hy)';     % y-grid
z  = ((1:modesizes(3))*hz - 0.5*hz)';     % z-grid
v1 = ((1:modesizes(4))*hv1 - 0.5*hv1 - 2*pi)';  % v1-grid
v2 = ((1:modesizes(5))*hv2 - 0.5*hv2 - 2*pi)';  % v2-grid
v3 = ((1:modesizes(6))*hv3 - 0.5*hv3 - 2*pi)';  % v3-grid

% === 3. Define BGK-type Maxwellian function ===
maxwellian_bgk = @(i,j,p,q,m,n) local_maxwellian(x(i), y(j), z(p), v1(q), v2(m), v3(n));

% === 4. Run HTACA compression ===
base_rel_tol = 1e-6;
tol_ratio = 0.1;
min_rank = ht.get_default_min_rank(tree);
max_rank = ht.get_default_max_rank(tree);
A = ht.HTACA(maxwellian_bgk, tree, base_rel_tol, tol_ratio, min_rank, max_rank);

disp('--- HTACA Compression Summary ---');
disp(['Final rank: ', mat2str(A.rank)]);
disp(['Frobenius norm: ', num2str(ht.norm(A))]);

% === 5. Evaluate accuracy on random indices ===
sample_size = 100;
sample_size = 100;
ind_mat = zeros(sample_size, 6);
for d = 1:6
    ind_mat(:, d) = randi(modesizes(d), sample_size, 1);
end
exact_vals = zeros(sample_size, 1);
approx_vals = zeros(sample_size, 1);

for k = 1:sample_size
    i = ind_mat(k, 1); j = ind_mat(k, 2); p = ind_mat(k, 3);
    q = ind_mat(k, 4); m = ind_mat(k, 5); n = ind_mat(k, 6);
    exact_vals(k)  = maxwellian_bgk(i,j,p,q,m,n);
    approx_vals(k) = ht.evaluate_index(A, i,j,p,q,m,n);
end

rel_error = norm(exact_vals - approx_vals) / norm(exact_vals);
max_error = max(abs(exact_vals - approx_vals));

disp('--- Evaluation Error ---');
disp(['Relative error on 100 samples: ', num2str(rel_error)]);
disp(['Max absolute error: ', num2str(max_error)]);

% === Local Maxwellian function ===
function val = local_maxwellian(x, y, z, vx, vy, vz)
    rho = 1 + 0.2 * sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);  % density
    u1  = 0.1 * sin(pi*x);
    u2  = 0.1 * sin(pi*y);
    u3  = 0.1 * sin(pi*z);
    T   = 1;
    v_rel_sq = (vx - u1).^2 + (vy - u2).^2 + (vz - u3).^2;
    val = rho ./ ((2*pi*T)^(3/2)) .* exp(-v_rel_sq / (2*T));
end
