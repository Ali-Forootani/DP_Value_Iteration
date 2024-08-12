function adjust_lambda_mu()
    % Parameters
    m = 4;  % 4 dimensions for lambda and mu
    N = 3;
    tolerance = 0.1;
    max_iterations = 100;  % Maximum number of iterations to avoid infinite loop

    % Generate initial random values for lambda and mu
    [lambda, mu] = generate_random_lambda_mu(m);

    % Construct the state space
    [S, ~] = construct_state_space(m, N);

    % Initialize variables to store best results
    best_lambda = lambda;
    best_mu = mu;
    best_max_diff = inf;
    best_state = [];

    % Loop to find suitable lambda and mu values
    for iteration = 1:max_iterations
        % Compute the maximum difference in transition probabilities
        [max_diff, state_with_max_diff] = max_probability_difference_all(S, N, m, lambda, mu);
        
        % If this iteration yields a lower max difference, update the best values
        if max_diff < best_max_diff
            best_max_diff = max_diff;
            best_lambda = lambda;
            best_mu = mu;
            best_state = state_with_max_diff;
        end
        
        % Check if the difference is within the tolerance
        if best_max_diff < tolerance
            fprintf('Suitable lambda and mu found:\n');
            fprintf('Lambda: [%f, %f, %f, %f]\n', best_lambda);
            fprintf('Mu: [%f, %f, %f, %f]\n', best_mu);
            fprintf('State with maximum difference: [%d, %d, %d, %d]\n', best_state);
            fprintf('Maximum difference in transition probabilities: %.4f\n', best_max_diff);
            save('best_params.mat', 'best_lambda', 'best_mu', 'best_max_diff', 'best_state');
            return;
        end
        
        % Generate new random lambda and mu values
        [lambda, mu] = generate_random_lambda_mu(m);
    end
    
    fprintf('Could not find suitable lambda and mu within %d iterations.\n', max_iterations);
    fprintf('Best Lambda: [%f, %f, %f, %f]\n', best_lambda);
    fprintf('Best Mu: [%f, %f, %f, %f]\n', best_mu);
    fprintf('State with maximum difference: [%d, %d, %d, %d]\n', best_state);
    fprintf('Best Maximum difference in transition probabilities: %.4f\n', best_max_diff);
    save('best_params.mat', 'best_lambda', 'best_mu', 'best_max_diff', 'best_state');
end

function [lambda, mu] = generate_random_lambda_mu(m)
    % Generate random values for lambda and mu such that lambda[i] + mu[i] < 1
    lambda = rand(m, 1) * 0.5;  % Start with random values between 0 and 0.5
    mu = rand(m, 1) * 0.5;      % Start with random values between 0 and 0.5

    for i = 1:m
        % Ensure that lambda[i] + mu[i] < 1
        while lambda(i) + mu(i) >= 1
            lambda(i) = rand * 0.5;
            mu(i) = rand * 0.5;
        end
    end
end

function [S, NS] = construct_state_space(m, N)
    % Construct the state space matrix S
    n = 0;
    l = 0;
    while n <= ((N+1)^m - 1)
        s = base(m, N, n);
        if sum(s) <= N
            l = l + 1;
            S(l, 1:m) = s;
        end
        n = n + 1;
    end
    NS = size(S);
end

function [max_diff, state_with_max_diff] = max_probability_difference_all(S, N, m, lambda, mu)
    % Compute the maximum difference in transition probabilities
    max_diff = 0;
    NS = size(S, 1);
    state_with_max_diff = [];

    for n1 = 1:NS
        for a1 = 1:m
            for a2 = a1+1:m
                if a1 <= size(lambda, 1) && a2 <= size(lambda, 1)
                    try
                        [transition1, combined_probability1] = stateanalysis(S(n1, :), N, m, a1, lambda, mu);
                        [transition2, combined_probability2] = stateanalysis(S(n1, :), N, m, a2, lambda, mu);

                        probabilities1 = compute_probabilities(S, transition1, combined_probability1);
                        probabilities2 = compute_probabilities(S, transition2, combined_probability2);

                        sum_prob1 = sum(probabilities1);
                        sum_prob2 = sum(probabilities2);

                        if abs(sum_prob1 - 1) > 1e-6
                            if sum_prob1 > 0
                                probabilities1 = probabilities1 / sum_prob1;
                            else
                                error(['Sum of probabilities1 for state ', num2str(n1), ' and action ', num2str(a1), ' is zero.']);
                            end
                        end

                        if abs(sum_prob2 - 1) > 1e-6
                            if sum_prob2 > 0
                                probabilities2 = probabilities2 / sum_prob2;
                            else
                                error(['Sum of probabilities2 for state ', num2str(n1), ' and action ', num2str(a2), ' is zero.']);
                            end
                        end

                        differences = abs(probabilities1 - probabilities2);
                        max_difference = max(differences);

                        if max_difference > max_diff
                            max_diff = max_difference;
                            state_with_max_diff = S(n1, :);
                        end

                    catch ME
                        disp(['Error processing state ', num2str(n1), ' with actions ', num2str(a1), ' and ', num2str(a2)]);
                        disp(ME.message);
                    end
                end
            end
        end
    end

    disp(['The maximum difference in transition probabilities: ', num2str(max_diff)]);
end

function probabilities = compute_probabilities(S, transition, combined_probability)
    % Initialize probabilities array
    NS = size(S, 1);
    probabilities = zeros(1, NS);

    % Fill the probabilities for the transitions
    for n2 = 1:size(transition, 1)
        idx = find(ismember(S, transition(n2, :), 'rows'));
        if ~isempty(idx)
            probabilities(idx) = prod(combined_probability(n2, :));
        end
    end
end

function s = base(m, N, n)
    % Generate state in base (N+1) representation
    s = zeros(1, m);
    for i = m:-1:1
        s(i) = mod(n, N + 1);
        n = floor(n / (N + 1));
    end
end
