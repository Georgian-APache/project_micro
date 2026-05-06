% =========================================================================
% elevation_sweep.m — Performance gain vs satellite elevation angle
% =========================================================================
% Sweeps elevation from 5° to 90° at a fixed SNR and shows how the RIS
% advantage grows at low elevation.
% =========================================================================

clear; clc; close all; rng(11);

% Parameters
N               = 64;
SNR_dB_fixed    = 5;            % fixed SNR for the sweep
elev_range_deg  = 5:5:85;       % elevation sweep
num_monte_carlo = 200;
num_bits        = 1000;

K_min_dB = -5;
K_max_dB = 12;
K_dB_ris = 10;

pso_swarm = 20; pso_iter = 30;

ber_no  = zeros(size(elev_range_deg));
ber_pso = zeros(size(elev_range_deg));

snr_lin   = 10^(SNR_dB_fixed/10);
noise_var = 1/snr_lin;

fprintf('Sweeping elevation at SNR = %d dB...\n', SNR_dB_fixed);
for e_idx = 1:length(elev_range_deg)
    elev = elev_range_deg(e_idx);
    K_dB_d = K_min_dB + (K_max_dB - K_min_dB) * sind(elev);

    err_no = 0; err_pso = 0; total = 0;
    for mc = 1:num_monte_carlo
        h_d  = rician_channel(1, K_dB_d);
        h_sr = rician_channel(N, K_dB_ris);
        h_rg = rician_channel(N, K_dB_ris);
        cascaded = h_sr .* h_rg;

        phi_pso = pso_optimize(h_d, cascaded, pso_swarm, pso_iter);
        h_eff_no  = h_d;
        h_eff_pso = h_d + sum(cascaded .* phi_pso);

        bits = randi([0 1], num_bits, 1);
        x    = 2*bits - 1;
        n    = sqrt(noise_var/2)*(randn(num_bits,1)+1j*randn(num_bits,1));

        b_no  = real((h_eff_no  * x + n) * conj(h_eff_no )) > 0;
        b_pso = real((h_eff_pso * x + n) * conj(h_eff_pso)) > 0;

        err_no  = err_no  + sum(b_no  ~= bits);
        err_pso = err_pso + sum(b_pso ~= bits);
        total   = total + num_bits;
    end
    ber_no(e_idx)  = err_no  / total;
    ber_pso(e_idx) = err_pso / total;
    fprintf('  elev=%2d°  NoRIS BER=%.4f  PSO BER=%.4f\n', ...
            elev, ber_no(e_idx), ber_pso(e_idx));
end

figure('Name', 'BER vs Elevation', 'Position', [100 100 800 500]);
semilogy(elev_range_deg, max(ber_no, 1e-6),  'r-o', 'LineWidth', 2); hold on;
semilogy(elev_range_deg, max(ber_pso,1e-6),  'b-^', 'LineWidth', 2);
grid on;
xlabel('Elevation Angle (degrees)', 'FontSize', 12);
ylabel('BER', 'FontSize', 12);
title(sprintf('BER vs Elevation Angle | SNR = %d dB | N = %d', ...
              SNR_dB_fixed, N), 'FontSize', 13);
legend({'No RIS', 'PSO-Optimized RIS'}, 'Location', 'NorthEast');
saveas(gcf, 'BER_vs_Elevation.png');
fprintf('Saved BER_vs_Elevation.png\n');
