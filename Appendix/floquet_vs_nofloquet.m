% Define the transition matrices for a PMDP with period τ=2
P0 = [0.8, 0.2; 0.1, 0.9];
P1 = [0.9, 0.1; 0.3, 0.7];

% Compute the monodromy matrix (product over one period)
P = P0 * P1;

% Compute the Floquet matrix as the square root of the monodromy matrix
tilde_P = sqrtm(P);

% Compute the transformation matrices Lambda_0 and Lambda_1
Lambda_0 = eye(size(P, 1));  % Identity matrix for Lambda_0
Lambda_1 = P0 * inv(tilde_P);  % Compute Lambda_1

% Define the reward vectors R0 and R1
R0 = [1; 0];  % Example reward vector for time step 0
R1 = [0; 1];  % Example reward vector for time step 1

% Transform the reward vectors
R0_transformed = inv(Lambda_0) * R0;
R1_transformed = inv(Lambda_1) * R1;

% Define discount factor
alpha = 0.9;

% Initialize value function at the final time step J_T (we assume it is known)
J_T = [0.0; 0.0];  % Example value function at the final time step

% Set the total number of time steps (we will calculate recursively from T to 0)
T = 50;  % Number of time steps to simulate
J_values_floquet = [];  % List to store value functions with Floquet transformation
J_values_no_floquet = [];  % List to store value functions without Floquet transformation

% Start the recursive Bellman equation from the final time step
J_k_floquet = J_T;
J_values_floquet = [J_values_floquet, J_k_floquet];  % Add the initial J_T to the list

J_k_no_floquet = J_T;
J_values_no_floquet = [J_values_no_floquet, J_k_no_floquet];  % Add the initial J_T to the list

% Recursively compute J_k for each time step with the Floquet transformation
for t = T-1:-1:0
    % With Floquet transformation (using tilde_P and transformed rewards)
    if mod(t, 2) == 0  % For even time steps, use R0_transformed
        J_k_floquet = R0_transformed + alpha * (tilde_P * J_k_floquet);
    else  % For odd time steps, use R1_transformed
        J_k_floquet = R1_transformed + alpha * (tilde_P * J_k_floquet);
    end
    
    J_values_floquet = [J_values_floquet, J_k_floquet];
    
    % Without Floquet transformation (using original transition matrices P0, P1)
    if mod(t, 2) == 0  % For even time steps, use R0
        J_k_no_floquet = R0 + alpha * (P0 * J_k_no_floquet);
    else  % For odd time steps, use R1
        J_k_no_floquet = R1 + alpha * (P1 * J_k_no_floquet);
    end
    
    J_values_no_floquet = [J_values_no_floquet, J_k_no_floquet];
end

% Reverse the order of the time steps for backward iteration
time_steps = T:-1:0;

% Plot the evolution of J_k over time for both methods (backward iteration)
figure('Position', [100, 100, 800, 550]); set(gcf, 'Color', 'w');

% Plot value functions for the state 1
plot(time_steps, J_values_floquet(1, :), '-', 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', [0.2, 0.4, 0.8], 'DisplayName', 'State 1 (Floquet)'); hold on;
plot(time_steps, J_values_no_floquet(1, :), '--', 'Color', [0.1, 0.5, 0.7], 'LineWidth', 2, 'Marker', 's', 'MarkerFaceColor', [0.1, 0.5, 0.7], 'DisplayName', 'State 1 (No Floquet)');

% Plot value functions for the state 2
plot(time_steps, J_values_floquet(2, :), '-', 'Color', [1, 0.4, 0.2], 'LineWidth', 2, 'Marker', '^', 'MarkerFaceColor', [1, 0.4, 0.2], 'DisplayName', 'State 2 (Floquet)');
plot(time_steps, J_values_no_floquet(2, :), '--', 'Color', [0.8, 0.1, 0.2], 'LineWidth', 2, 'Marker', 'd', 'MarkerFaceColor', [0.8, 0.1, 0.2], 'DisplayName', 'State 2 (No Floquet)');

% Adding labels and title
xlabel('Time Step (k)', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');


% --- PLOT: ΔJ_k(x) ---
% For simplicity, let's compute the residuals ΔJ for both methods
Delta_J_floquet = diff(J_values_floquet, 1, 2);  % Residuals for Floquet method
Delta_J_no_floquet = diff(J_values_no_floquet, 1, 2);  % Residuals for No Floquet method

figure('Position', [100, 100, 800, 500]); set(gcf, 'Color', 'w');
% Plot value functions for the state 1
plot(time_steps, J_values_floquet(1, :), '-', 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', [0.2, 0.4, 0.8], 'DisplayName', 'State 1 (Floquet)'); hold on;
plot(time_steps, J_values_no_floquet(1, :), '--', 'Color', [0.1, 0.5, 0.7], 'LineWidth', 2, 'Marker', 's', 'MarkerFaceColor', [0.1, 0.5, 0.7], 'DisplayName', 'State 1 (No Floquet)');

% Plot value functions for the state 2
plot(time_steps, J_values_floquet(2, :), '-', 'Color', [1, 0.4, 0.2], 'LineWidth', 2, 'Marker', '^', 'MarkerFaceColor', [1, 0.4, 0.2], 'DisplayName', 'State 2 (Floquet)');
plot(time_steps, J_values_no_floquet(2, :), '--', 'Color', [0.8, 0.1, 0.2], 'LineWidth', 2, 'Marker', 'd', 'MarkerFaceColor', [0.8, 0.1, 0.2], 'DisplayName', 'State 2 (No Floquet)');

% Adding labels and title
xlabel('Time Step (k)', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');

