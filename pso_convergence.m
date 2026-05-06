% =========================================================================
% pso_convergence.m — Visualize PSO convergence behavior
% =========================================================================
% Generates a single random channel realization and plots how PSO's best
% fitness evolves over iterations. Compares against the closed-form optimum.
% =========================================================================

clear; clc; close all; rng(7);

% Parameters
N             = 64;
K_dB_direct   = -3;
K_dB_ris      = 10;
swarm_size    = 20;
max_iter      = 100;
num_trials    = 5;     % run multiple trials to show variability

% Generate ONE channel realization (shared across trials)
h_d  = rician_channel(1, K_dB_direct);
h_sr = rician_channel(N, K_dB_ris);
h_rg = rician_channel(N, K_dB_ris);
cascaded = h_sr .* h_rg;

% Closed-form optimum (theoretical ceiling)
phi_opt = exp(1j * (angle(h_d) - angle(cascaded)));
optimal_fitness = abs(h_d + sum(cascaded .* phi_opt))^2;

% Run PSO multiple times to show variability
all_histories = zeros(num_trials, max_iter);
for t = 1:num_trials
    [~, history] = pso_optimize(h_d, cascaded, swarm_size, max_iter);
    all_histories(t, :) = history;
end

% Plot
figure('Name', 'PSO Convergence', 'Position', [100 100 800 500]);
plot(1:max_iter, all_histories.', 'LineWidth', 1.2); hold on;
plot(1:max_iter, mean(all_histories,1), 'b-', 'LineWidth', 2.5);
yline(optimal_fitness, 'r--', 'LineWidth', 2, ...
      'Label', 'Closed-form optimum', 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('PSO Iteration', 'FontSize', 12);
ylabel('Best |h_{eff}|^2 (fitness)', 'FontSize', 12);
title(sprintf('PSO Convergence | N = %d Elements | %d Independent Runs', ...
              N, num_trials), 'FontSize', 13);
legend([repmat({'Trial'}, 1, num_trials), {'Mean across trials'}, ...
        {'Theoretical optimum'}], 'Location', 'SouthEast');

saveas(gcf, 'PSO_Convergence.png');
fprintf('Saved PSO_Convergence.png\n');
fprintf('Optimum fitness        : %.4f\n', optimal_fitness);
fprintf('Mean PSO final fitness : %.4f\n', mean(all_histories(:, end)));
fprintf('PSO achieved %.1f%% of optimum on average.\n', ...
        100 * mean(all_histories(:, end)) / optimal_fitness);
