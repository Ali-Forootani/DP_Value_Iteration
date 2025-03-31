% Apply consistent font settings to all axes
set(findall(gcf, 'Type', 'axes'), 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');

% Define the transition matrices
P0 = [0.8, 0.2; 0.1, 0.9];
P1 = [0.9, 0.1; 0.3, 0.7];

% Compute the monodromy matrix P (product over one period)
P_monodromy = P0 * P1;

% Compute the square root of the monodromy matrix (Floquet matrix)
P_sqrt = fractional_matrix_power(P_monodromy, 0.5);

% Compute the transformation matrix Lambda_1 using the formula
Lambda_0 = eye(2);  % Identity matrix
Lambda_1 = P0 / P_sqrt;  % Lambda_1 = P0 * inv(P_sqrt)

% Define the initial reward vectors
R0 = [1; 0];  % Reward vector at time step 0
R1 = [0; 1];  % Reward vector at time step 1

% Compute the transformed reward vectors
tilde_R0 = inv(Lambda_0) * R0;
tilde_R1 = inv(Lambda_1) * R1;

% Define the discount factor alpha
alpha = 0.9;

% Initialize the value function J as a 2x1 zero vector at the final step (T=50)
J = zeros(2, 1);

% Store the values of J for plotting
J_values = [];

% Number of iterations (steps) and starting at T=50
num_iterations = 50;

% Define the backward Bellman operator for a given time step
for t = num_iterations-1:-1:0  % Iterating backward
    % Determine the current period based on t (periodicity)
    period_index = mod(t, 2) + 1;
    if period_index == 1
        P_k = P0;
        R_k = R0;
    else
        P_k = P1;
        R_k = R1;
    end

    % Apply the backward Bellman operator for this period
    J = backward_bellman_operator(J, P_k, R_k, alpha);
    
    % Append the current value of J to the list
    J_values = [J_values; J'];
end

% Reverse the elements in J_values
J_values_reversed = flipud(J_values);

% Plot the behavior of J (both components of the value function) after reversing the elements
figure;
plot(J_values_reversed(:, 1), 'b-o', 'DisplayName', 'State 1');
hold on;
plot(J_values_reversed(:, 2), 'r-s', 'DisplayName', 'State 2');

% Apply consistent font settings to all axes
set(findall(gcf, 'Type', 'axes'), 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
% Highlight the first value with a different marker
%firstValue1 = J_values_reversed(1,1);
%firstValue2 = J_values_reversed(1,2);
%plot(1, firstValue1, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Highlight first value of J[0]
%plot(1, firstValue2, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Highlight first value of J[1]

hold off;

% Label settings
xlabel('Time Step', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Value Function', 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');

legend('show');
grid on;

% Fractional matrix power function
function result = fractional_matrix_power(A, power)
    [V, D] = eig(A);  % Eigenvalue decomposition
    D_new = D.^power;  % Raise the diagonal matrix to the desired power
    result = V * D_new / V;  % Reconstruct the matrix
end

% Backward Bellman operator function
function J_new = backward_bellman_operator(J, P, R, alpha)
    % Apply the backward Bellman operator to compute the value function for the given state.
    J_new = R + alpha * P' * J;  % Using P' for backward computation
end
