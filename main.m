% =========================================================================
% main.m — RIS-Assisted CubeSat-to-Ground Communication Simulation
% =========================================================================
% Compares four configurations:
%   1. No RIS                      (baseline)
%   2. Random RIS phases           (lower bound)
%   3. PSO-optimized RIS phases    (proposed algorithm)
%   4. Closed-form optimal phases  (upper bound)
%
% Output: BER vs SNR plot for all four configurations.
%
% Author: APache  | Project: RIS-Assisted CubeSat Communication
% =========================================================================

clear; clc; close all; rng(42);  % rng(42) for reproducibility

% -------------------------------------------------------------------------
% SIMULATION PARAMETERS
% -------------------------------------------------------------------------
N               = 64;             % Number of RIS elements
elevation_deg   = 10;             % Low elevation angle (degrees)
SNR_dB_range    = -10:2:20;       % SNR sweep (dB)
num_monte_carlo = 200;            % Monte Carlo iterations per SNR
num_bits        = 1000;           % Bits per Monte Carlo iteration

% Elevation-dependent Rician K-factors (in dB)
K_min_dB = -5;   % at 0 deg elevation (worst case)
K_max_dB = 12;   % at 90 deg elevation (best case)
K_dB_direct = K_min_dB + (K_max_dB - K_min_dB) * sind(elevation_deg);
K_dB_ris    = 10;   % RIS sub-links assumed strong LOS (mounted high, clear view)

fprintf('Direct path K-factor at %d deg elevation: %.2f dB\n', ...
        elevation_deg, K_dB_direct);

% PSO parameters
pso_swarm_size  = 20;
pso_iterations  = 30;

% -------------------------------------------------------------------------
% RESULTS STORAGE
% -------------------------------------------------------------------------
ber_no_ris      = zeros(size(SNR_dB_range));
ber_random_ris  = zeros(size(SNR_dB_range));
ber_pso_ris     = zeros(size(SNR_dB_range));
ber_optimal_ris = zeros(size(SNR_dB_range));

% -------------------------------------------------------------------------
% MAIN MONTE CARLO LOOP
% -------------------------------------------------------------------------
fprintf('\nRunning simulation...\n');
tic;
for snr_idx = 1:length(SNR_dB_range)
    snr_db     = SNR_dB_range(snr_idx);
    snr_linear = 10^(snr_db/10);
    noise_var  = 1/snr_linear;

    err_no = 0; err_rand = 0; err_pso = 0; err_opt = 0;
    total_bits = 0;

    for mc = 1:num_monte_carlo
        % --- Generate channel realizations ---
        h_d  = rician_channel(1, K_dB_direct);
        h_sr = rician_channel(N, K_dB_ris);
        h_rg = rician_channel(N, K_dB_ris);

        % Cascaded channel (per-element product)
        cascaded = h_sr .* h_rg;     % N x 1

        % --- Phase configurations ---
        % 1. No RIS:    phi = 0 (RIS not present)
        % 2. Random:    uniform in [0, 2*pi)
        % 3. Optimal:   close-form alignment with h_d
        % 4. PSO:       iterative optimizer
        phi_random  = exp(1j * 2*pi*rand(N,1));
        phi_optimal = exp(1j * (angle(h_d) - angle(cascaded)));
        phi_pso     = pso_optimize(h_d, cascaded, ...
                                   pso_swarm_size, pso_iterations);

        % --- Effective channels ---
        h_eff_no   = h_d;
        h_eff_rand = h_d + sum(cascaded .* phi_random);
        h_eff_pso  = h_d + sum(cascaded .* phi_pso);
        h_eff_opt  = h_d + sum(cascaded .* phi_optimal);

        % --- BPSK transmission ---
        bits  = randi([0 1], num_bits, 1);
        x     = 2*bits - 1;                              % +1 / -1
        noise = sqrt(noise_var/2) * (randn(num_bits,1) + 1j*randn(num_bits,1));

        y_no   = h_eff_no   * x + noise;
        y_rand = h_eff_rand * x + noise;
        y_pso  = h_eff_pso  * x + noise;
        y_opt  = h_eff_opt  * x + noise;

        % --- Coherent BPSK demodulation ---
        bits_no   = real(y_no   * conj(h_eff_no  )) > 0;
        bits_rand = real(y_rand * conj(h_eff_rand)) > 0;
        bits_pso  = real(y_pso  * conj(h_eff_pso )) > 0;
        bits_opt  = real(y_opt  * conj(h_eff_opt )) > 0;

        % --- Count errors ---
        err_no   = err_no   + sum(bits_no   ~= bits);
        err_rand = err_rand + sum(bits_rand ~= bits);
        err_pso  = err_pso  + sum(bits_pso  ~= bits);
        err_opt  = err_opt  + sum(bits_opt  ~= bits);
        total_bits = total_bits + num_bits;
    end

    ber_no_ris(snr_idx)      = err_no   / total_bits;
    ber_random_ris(snr_idx)  = err_rand / total_bits;
    ber_pso_ris(snr_idx)     = err_pso  / total_bits;
    ber_optimal_ris(snr_idx) = err_opt  / total_bits;

    fprintf('  SNR = %3d dB | NoRIS=%.4f | Rand=%.4f | PSO=%.4f | Opt=%.4f\n', ...
            snr_db, ber_no_ris(snr_idx), ber_random_ris(snr_idx), ...
            ber_pso_ris(snr_idx), ber_optimal_ris(snr_idx));
end
elapsed = toc;
fprintf('\nSimulation complete in %.1f seconds.\n', elapsed);

% -------------------------------------------------------------------------
% PLOT BER vs SNR
% -------------------------------------------------------------------------
figure('Name', 'BER vs SNR', 'Position', [100 100 800 600]);
semilogy(SNR_dB_range, max(ber_no_ris,      1e-6), 'r-o',  'LineWidth', 2); hold on;
semilogy(SNR_dB_range, max(ber_random_ris,  1e-6), 'k--s', 'LineWidth', 2);
semilogy(SNR_dB_range, max(ber_pso_ris,     1e-6), 'b-^',  'LineWidth', 2);
semilogy(SNR_dB_range, max(ber_optimal_ris, 1e-6), 'g-d',  'LineWidth', 2);
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title(sprintf('BER vs SNR | %d° Elevation | N = %d RIS Elements', ...
              elevation_deg, N), 'FontSize', 13);
legend({'No RIS (baseline)', 'Random RIS (lower bound)', ...
        'PSO-optimized RIS', 'Closed-form optimal (upper bound)'}, ...
       'Location', 'SouthWest', 'FontSize', 11);
ylim([1e-5 1]);

% -------------------------------------------------------------------------
% PRINT QUANTITATIVE GAIN AT BER = 1e-3
% -------------------------------------------------------------------------
target_ber = 1e-3;
fprintf('\n=== Performance Summary at BER = %g ===\n', target_ber);
snr_no  = interp_snr(SNR_dB_range, ber_no_ris,      target_ber);
snr_pso = interp_snr(SNR_dB_range, ber_pso_ris,     target_ber);
snr_opt = interp_snr(SNR_dB_range, ber_optimal_ris, target_ber);
fprintf('  No RIS         : SNR required ≈ %.1f dB\n', snr_no);
fprintf('  PSO-RIS        : SNR required ≈ %.1f dB\n', snr_pso);
fprintf('  Optimal RIS    : SNR required ≈ %.1f dB\n', snr_opt);
fprintf('  PSO-RIS gain over No-RIS  : %.1f dB\n', snr_no - snr_pso);
fprintf('  PSO-to-Optimal gap        : %.1f dB\n', snr_pso - snr_opt);

% Save figure
saveas(gcf, 'BER_vs_SNR.png');
fprintf('\nFigure saved as BER_vs_SNR.png\n');

% -------------------------------------------------------------------------
% HELPER: interpolate SNR required to reach a target BER
% -------------------------------------------------------------------------
function snr_req = interp_snr(snr_vec, ber_vec, target)
    valid = ber_vec > 0 & ber_vec < 1;
    if nnz(valid) < 2
        snr_req = NaN; return;
    end
    log_ber = log10(ber_vec(valid));
    snr_v   = snr_vec(valid);
    if target < min(ber_vec(valid)) || target > max(ber_vec(valid))
        snr_req = NaN; return;
    end
    snr_req = interp1(log_ber, snr_v, log10(target), 'linear');
end
