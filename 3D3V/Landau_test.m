% === Add necessary folders ===
addpath('Methed Related');
addpath('ht_toolbox');
addpath('SetUp Initialization');
addpath('Gemini');

for i = 1:769
    rho_ht = compute_rho(f_ht{i},v_modes,h_mesh);
    E_field = get_E_field_FFT(rho_ht,dom_len,h_mesh,base_rel_eps,tol_ratio,min_rank,max_rank,kx);
    E_amplitude(i) = real(compute_electric_field_energy(h_mesh(v_modes), E_field));
end
% t = t(1:155);
% E_amplitude = abs(normphi);
log_E = log(E_amplitude + 1e-10);
% Load the saved electric field data

% Select a time range where damping occurs (e.g., after initial transient)
t_start_idx = 15;  % Skip initial transients
t_end_idx = floor(length(t)/2);  % Use half of the time series

exclude_indices = [];
t_fit = t(t_start_idx:t_end_idx);
log_E_fit = log_E(t_start_idx:t_end_idx);
t_fit(exclude_indices) = [];
log_E_fit(exclude_indices) = [];

% Perform a linear fit to extract damping rate
coeffs = polyfit(t_fit, log_E_fit, 1);
gamma = coeffs(1);  % Landau damping rate (negative slope)

% Generate fitted curve
log_E_fit_line = polyval(coeffs, t)+0.8;

% Plot Landau damping
figure;
plot(t, log_E, 'b','LineWidth', 2, 'DisplayName', 'Log |E|');
hold on;
plot(t, log_E_fit_line, 'r--','LineWidth', 2, 'DisplayName', sprintf('Fit: \\gamma = %.4f', gamma));

xlabel('Time', 'FontSize', 16);
ylabel('electrostatic energy', 'FontSize', 16);
% title('Landau Damping Rate Extraction', 'FontSize', 18);
legend('FontSize', 14);
grid on;
set(gca, 'FontSize', 14);  % Set font size for axes
hold off;

fprintf('Estimated Landau damping rate γ = %.4f\n', gamma);

