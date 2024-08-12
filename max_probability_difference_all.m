function max_diff = max_probability_difference_all(S, N, m, lambda, mu)
    % Initialize the maximum difference variable
    max_diff = 0;

    % Get the number of states in S
    NS = size(S, 1);

    % Loop through each state in S
    for n1 = 1:NS
        % For each state, compare probabilities for every pair of actions
        for a1 = 1:m
            for a2 = a1+1:m
                % Ensure action indices are within bounds
                if a1 <= size(lambda, 1) && a2 <= size(lambda, 1)
                    % Compute transition probabilities for action a1
                    try
                        [transition1, combined_probability1] = stateanalysis(S(n1, :), N, m, a1, lambda, mu);
                        [transition2, combined_probability2] = stateanalysis(S(n1, :), N, m, a2, lambda, mu);

                        % Initialize probability arrays for entire state space
                        probabilities1 = zeros(1, NS);
                        probabilities2 = zeros(1, NS);

                        % Fill the probabilities for transitions under action a1
                        for n2 = 1:size(transition1, 1)
                            % Find the index of the corresponding state in S
                            idx = find(ismember(S, transition1(n2, :), 'rows'));
                            probabilities1(idx) = prod(combined_probability1(n2, :));
                        end

                        % Fill the probabilities for transitions under action a2
                        for n2 = 1:size(transition2, 1)
                            % Find the index of the corresponding state in S
                            idx = find(ismember(S, transition2(n2, :), 'rows'));
                            probabilities2(idx) = prod(combined_probability2(n2, :));
                        end
                        
                       

                        % Check if the sum of probabilities is 1, normalize if necessary
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



                        % Compute the differences between the probabilities
                        differences = abs(probabilities1 - probabilities2);

                        % Update the maximum difference
                        max_diff = max(max_diff, max(differences));

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
