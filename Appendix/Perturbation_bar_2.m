clear; clc;

% MDP parameters
states = [0, 1];
actions = {'a', 'b'};
alpha = 0.9;
epsilon = 0.12;
beta = 0.2;
tau = 80;
H = 500;

% Rewards
R = containers.Map({0, 1}, [0.3, 0.1]);

% Transition probabilities
P_bar = containers.Map;
P_bar('0a') = [0.7, 0.3];
P_bar('1a') = [0.2, 0.8];
P_bar('0b') = [0.6, 0.4];
P_bar('1b') = [0.5, 0.5];

P_0 = containers.Map;
P_0('0a') = [0.6, 0.4];
P_0('1a') = [0.3, 0.7];
P_0('0b') = [0.5, 0.5];
P_0('1b') = [0.4, 0.6];

% Compute value functions
J_backward = backward_bellman_iteration(H, true, P_bar, P_0, R, alpha, epsilon, tau);
J_stationary = backward_bellman_iteration(H, false, P_bar, P_0, R, alpha, epsilon, tau);
J_bar = J_stationary(1, :);

% Compute ΔJ_k
Delta_J = J_backward - J_bar;

% Steady-state window
steady_start = 100;
steady_end = 200;
Delta_J_steady = Delta_J(steady_start:steady_end, :);

% Empirical max deviation
max_delta_steady = max(abs(Delta_J_steady), [], 1);
fprintf('Empirical max |ΔJ_k(x)| (steady-state k=100 to 200):\n');
for x = 1:2
    fprintf('  State %d: %.4f\n', x-1, max_delta_steady(x));
end

% Conservative theoretical bound
J_bar_inf = max(J_bar);
term1 = alpha * epsilon * beta * J_bar_inf;
fprintf('\nConservative base bound term: %.4f\n', term1);

% FFT and harmonic term
A1_values = zeros(1, 2);
for x = 1:2
    fft_vals = fft(Delta_J_steady(:, x));
    A1 = 2 * abs(fft_vals(2)) / size(Delta_J_steady, 1);
    A1_values(x) = A1;
    fprintf('A₁ estimate (state %d): %.4f\n', x-1, A1);
end

term2 = sqrt(2) * alpha * (1 + epsilon * beta) * max(A1_values);
fprintf('First harmonic term: %.4f\n', term2);

total_bound = term1 + term2;
fprintf('\nTotal theoretical bound: %.4f\n', total_bound);

% ==== PLOT ====
figure('Position', [100, 100, 700, 500]);
set(gcf, 'Color', 'w');

% Custom colors
color_J0 = [0 0 0.5];         % navy blue
color_J1 = [0.6 0 0];         % dark red
color_steady = [0 0.5 0];     % dark green
color_upper = [1 0.5 0];      % orange
color_lower = [0.5 0 0.5];    % purple

% Plot value functions
plot(J_backward(:, 1), '-', 'Color', color_J0, 'LineWidth', 2, 'DisplayName', 'J_k(0)');
hold on;
plot(J_backward(:, 2), '--', 'Color', color_J1, 'LineWidth', 2, 'DisplayName', 'J_k(1)');

% Plot steady-state values using overbar Unicode character
yline(J_bar(1), ':', 'Color', color_steady, 'LineWidth', 2, 'DisplayName', 'Steady-State J̄(0)');
yline(J_bar(2), ':', 'Color', color_steady, 'LineWidth', 2, 'HandleVisibility','off');

% Plot bounds
yline(J_bar(1) + total_bound, '--', 'Color', color_upper, 'LineWidth', 2, 'DisplayName', 'Upper Bound');
yline(J_bar(1) - total_bound, '--', 'Color', color_lower, 'LineWidth', 2, 'DisplayName', 'Lower Bound');
yline(J_bar(2) + total_bound, '--', 'Color', color_upper, 'LineWidth', 2, 'HandleVisibility','off');
yline(J_bar(2) - total_bound, '--', 'Color', color_lower, 'LineWidth', 2, 'HandleVisibility','off');

% Axis labels
xlabel('Time Step (k)', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');

legend('Location', 'best');
grid on;

% Apply font settings globally
set(findall(gcf, 'Type', 'axes'), 'FontName', 'Times New Roman', ...
    'FontSize', 20, 'FontWeight', 'bold');

% ======= FUNCTION DEFINITIONS BELOW =======

function P_k = compute_perturbed_transitions(k, epsilon, P_bar, P_0, tau)
    keys = P_bar.keys;
    P_k = containers.Map;
    f_val = abs(epsilon * sin(2 * pi * k / tau));
    for i = 1:length(keys)
        key = keys{i};
        P_k(key) = P_bar(key) + f_val * (P_0(key) - P_bar(key));
    end
end

function J = backward_bellman_iteration(H, use_perturbation, P_bar, P_0, R, ...
                                        alpha, epsilon, tau)
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
                if u{1} == 'a'
                    key = sprintf('%da', x);
                else
                    key = sprintf('%db', x);
                end
                val = max(val, R(x) + alpha * dot(P_k(key), J(k+1, :)));
            end
            J(k, x+1) = val;
        end
    end
end
