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
A_0 = 0;
A1_values = zeros(1, 2);
B1_values = zeros(1, 2);
N = size(Delta_J_steady, 1);
for x = 1:2
    fft_vals = fft(Delta_J_steady(:, x));
    a1 = 2 * real(fft_vals(2)) / N;
    b1 = 2 * imag(fft_vals(2)) / N;
    A1_values(x) = abs(a1);
    B1_values(x) = abs(b1);
end
A1_max = max(A1_values);
B1_max = max(B1_values);

bound = alpha * (1 + epsilon * beta) * (A_0 + A1_max + B1_max) + epsilon * beta * J_bar_inf;
fprintf('\nTotal theoretical bound: %.4f\n', bound);

% ==== PLOT ====
figure('Position', [100, 100, 700, 500]);
set(gcf, 'Color', 'w');

% Color settings
color_J0 = [0.1 0.2 0.7];       % Blue
color_J1 = [0.8 0.1 0.1];       % Red
color_bar = [0.1 0.6 0.1];      % Green
color_upper = [0.85 0.5 0.1];   % Orange
color_lower = [0.5 0.1 0.6];    % Purple

% Plot value functions
plot(J_backward(:, 1), '-', 'Color', color_J0, 'LineWidth', 2, 'DisplayName', 'J_k(0)');
hold on;
plot(J_backward(:, 2), '-', 'Color', color_J1, 'LineWidth', 2, 'DisplayName', 'J_k(1)');

% Plot steady-state value
yline(J_bar(1), ':', 'Color', color_bar, 'LineWidth', 2, 'DisplayName', 'Steady-State J̄(0)');
yline(J_bar(2), ':', 'Color', color_bar, 'LineWidth', 2, 'HandleVisibility', 'off');

% Plot bounds
yline(J_bar(1) + bound, '--', 'Color', color_upper, 'LineWidth', 2, 'DisplayName', 'Upper Bound');
yline(J_bar(1) - bound, '--', 'Color', color_lower, 'LineWidth', 2, 'DisplayName', 'Lower Bound');
yline(J_bar(2) + bound, '--', 'Color', color_upper, 'LineWidth', 2, 'HandleVisibility','off');
yline(J_bar(2) - bound, '--', 'Color', color_lower, 'LineWidth', 2, 'HandleVisibility','off');

xlabel('Time Step (k)', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');

legend('Location', 'best');
grid on;
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
