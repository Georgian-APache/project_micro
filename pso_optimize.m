function [phi_best, fitness_history] = pso_optimize(h_d, cascaded, swarm_size, max_iter)
% PSO_OPTIMIZE  Particle Swarm Optimization for RIS phase shifts.
%
% Maximizes |h_d + sum_n cascaded_n * exp(j*theta_n)|^2
% over phase vector theta in [0, 2*pi)^N.
%
% Inputs:
%   h_d         — direct channel coefficient (scalar complex)
%   cascaded    — N x 1 cascaded channel (h_sr .* h_rg)
%   swarm_size  — number of particles
%   max_iter    — number of PSO iterations
%
% Outputs:
%   phi_best        — N x 1 vector of optimal reflection coefficients (e^j*theta)
%   fitness_history — 1 x max_iter vector of best fitness per iteration

    N = length(cascaded);

    % PSO hyperparameters
    w  = 0.7;        % inertia weight
    c1 = 1.5;        % cognitive coefficient
    c2 = 1.5;        % social coefficient
    v_max = pi/4;    % max velocity per dimension

    % Initialize positions (phase angles) and velocities
    positions  = 2*pi*rand(swarm_size, N);
    velocities = zeros(swarm_size, N);

    % Evaluate initial fitness
    fitness = zeros(swarm_size, 1);
    for s = 1:swarm_size
        phi = exp(1j * positions(s,:).');
        fitness(s) = abs(h_d + sum(cascaded .* phi))^2;
    end

    pbest_pos = positions;
    pbest_fit = fitness;
    [gbest_fit, idx] = max(fitness);
    gbest_pos = positions(idx, :);

    fitness_history = zeros(1, max_iter);

    % Main PSO loop
    for iter = 1:max_iter
        for s = 1:swarm_size
            r1 = rand(1, N);
            r2 = rand(1, N);
            velocities(s,:) = w*velocities(s,:) ...
                + c1*r1.*(pbest_pos(s,:) - positions(s,:)) ...
                + c2*r2.*(gbest_pos      - positions(s,:));
            % Clamp velocity
            velocities(s,:) = max(min(velocities(s,:), v_max), -v_max);
            % Update position (wrap to [0, 2*pi))
            positions(s,:) = mod(positions(s,:) + velocities(s,:), 2*pi);

            % Evaluate fitness
            phi = exp(1j * positions(s,:).');
            fit = abs(h_d + sum(cascaded .* phi))^2;

            % Update personal best
            if fit > pbest_fit(s)
                pbest_fit(s)    = fit;
                pbest_pos(s,:)  = positions(s,:);
                % Update global best
                if fit > gbest_fit
                    gbest_fit = fit;
                    gbest_pos = positions(s,:);
                end
            end
        end
        fitness_history(iter) = gbest_fit;
    end

    phi_best = exp(1j * gbest_pos.');
end
