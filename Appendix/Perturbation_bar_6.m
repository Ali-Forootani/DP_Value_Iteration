clear; clc;

% ====================================
% Residual Value Analysis (Fourier Perturbation)
% ====================================

% Parameters
states = [0, 1];
actions = {'a', 'b'};
alpha = 0.9;
epsilon = 0.02;
tau = 80;
H = 500;

% Rewards
R = containers.Map({0, 1}, [0.3, 0.1]);

% Nominal transition probabilities
P_bar = containers.Map;
P_bar('0a') = [0.7, 0.3];
P_bar('1a') = [0.2, 0.8];
P_bar('0b') = [0.6, 0.4];
P_bar('1b') = [0.5, 0.5];

% Perturbed transition probabilities
P_0 = containers.Map;
P_0('0a') = [0.6, 0.4];
P_0('1a') = [0.3, 0.7];
P_0('0b') = [0.5, 0.5];
P_0('1b') = [0.4, 0.6];

% Compute value functions
J_full = backward_bellman_iteration(H, true, P_bar, P_0, R, alpha, epsilon, tau);
J_nominal = backward_bellman_iteration(H, false, P_bar, P_0, R, alpha, epsilon, tau);
J_bar = J_nominal(1, :);
Delta_J = J_full - J_bar;

% Steady-state window
steady_start = 100;
steady_end = 200;
Delta_J_steady = Delta_J(steady_start:steady_end, :);
max_delta_steady = max(abs(Delta_J_steady), [], 1);

fprintf('Empirical max |ΔJ_k(x)| (steady-state k=100 to 200):\n');
for x = 1:2
    fprintf('  State %d: %.4f\n', x-1, max_delta_steady(x));
end

% η(x) computation
eta = zeros(1, 2);
for x = 0:1
    delta_P = abs(P_0(sprintf('%da', x)) - P_bar(sprintf('%da', x)));
    eta(x+1) = dot(delta_P, J_bar);
end
J_max = max(abs(J_bar));
beta = max([sum(abs(P_0('0a') - P_bar('0a'))), sum(abs(P_0('1a') - P_bar('1a')))]);
bound = (alpha * epsilon * beta * J_max) / (1 - alpha);

fprintf('\nMax |η(x)| = %.6f\n', max(abs(eta)));
fprintf('Theoretical bound = %.6f\n\n', bound);

fprintf('η(x) values:\n');
for x = 0:1
    fprintf('  η(%d) = %.6f\n', x, eta(x+1));
end

% Analytical Delta J
k_vals = 0:H-1;
Delta_J_analytical = zeros(H, 2);
for x = 1:2
    Delta_J_analytical(:, x) = alpha * eta(x) * f_k(k_vals, epsilon, tau);
end

% ==== PLOT: J_k(x) with Bound ====
figure('Position', [100, 100, 800, 550]); set(gcf, 'Color', 'w');
color_J0 = [0.1 0.2 0.7]; color_J1 = [0.8 0.1 0.1];
color_bar = [0.1 0.6 0.1]; color_upper = [0.85 0.5 0.1]; color_lower = [0.5 0.1 0.6];

plot(J_full(:, 1), '-', 'Color', color_J0, 'LineWidth', 2, 'DisplayName', 'J_k(0)'); hold on;
plot(J_full(:, 2), '-', 'Color', color_J1, 'LineWidth', 2, 'DisplayName', 'J_k(1)');
yline(J_bar(1), ':', 'Color', color_bar, 'LineWidth', 2, 'DisplayName', 'Steady-State J̄(0)');
yline(J_bar(2), ':', 'Color', color_bar, 'LineWidth', 2, 'HandleVisibility', 'off');
yline(J_bar(1) + bound, '--', 'Color', color_upper, 'LineWidth', 2, 'DisplayName', 'Upper Bound');
yline(J_bar(1) - bound, '--', 'Color', color_lower, 'LineWidth', 2, 'DisplayName', 'Lower Bound');
yline(J_bar(2) + bound, '--', 'Color', color_upper, 'LineWidth', 2, 'HandleVisibility', 'off');
yline(J_bar(2) - bound, '--', 'Color', color_lower, 'LineWidth', 2, 'HandleVisibility', 'off');

xlabel('Time Step (k)', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');

% ==== PLOT: ΔJ_k(x) ====
figure('Position', [100, 100, 800, 500]); set(gcf, 'Color', 'w');
plot(Delta_J(:, 1), '-', 'Color', color_J0, 'LineWidth', 2, 'DisplayName', 'ΔJ_k(0) Sim'); hold on;
plot(Delta_J(:, 2), '-', 'Color', color_J1, 'LineWidth', 2, 'DisplayName', 'ΔJ_k(1) Sim');
plot(Delta_J_analytical(:, 1), '--', 'Color', color_J0, 'LineWidth', 2, 'DisplayName', 'ΔJ_k(0) Analytical');
plot(Delta_J_analytical(:, 2), '--', 'Color', color_J1, 'LineWidth', 2, 'DisplayName', 'ΔJ_k(1) Analytical');

xlabel('Time Step (k)', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Residual ΔJ_k(x)', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');

% ======= FUNCTION DEFINITIONS =======

function val = f_k(k, epsilon, tau)
    val = epsilon * sin(2 * pi * k / tau);
end

function P_k = compute_perturbed_transitions(k, epsilon, P_bar, P_0, tau)
    keys = P_bar.keys;
    P_k = containers.Map;
    f_val = f_k(k, epsilon, tau);
    for i = 1:length(keys)
        key = keys{i};
        P_k(key) = P_bar(key) + f_val * (P_0(key) - P_bar(key));
    end
end

function J = backward_bellman_iteration(H, use_perturbation, P_bar, P_0, R, alpha, epsilon, tau)
    J = zeros(H + 1, 2);
    for k = H:-1:1
        if use_perturbation
            P_k = compute_perturbed_transitions(mod(k, tau), epsilon, P_bar, P_0, tau);
        else
            P_k = P_bar;
        end
        for x = 0:1
            val = -inf;
            for u = {'a', 'b'}
                key = sprintf('%d%s', x, u{1});
                val = max(val, R(x) + alpha * dot(P_k(key), J(k+1, :)));
            end
            J(k, x+1) = val;
        end
    end
end
