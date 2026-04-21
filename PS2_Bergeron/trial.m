
% do iteration
while error > Params.e_stop * (1-Params.beta) % continue until convergence
    for i in 1:size(state_space, 1) % loop over all states
        a = state_space(i, 1); % current wealth
        z = state_space(i, 2); % current income

        
        % compute expected value for each possible action (consumption choice)
        expected_values = zeros(size(a_grid)); % initialize expected values for each action
        for j = 1:length(a_grid) % loop over possible actions (consumption choices)
            a_next = a_grid(j); % next period's wealth based on current action
            c = max(0, a + z - a_next); % consumption based on current wealth, income, and next period's wealth
            u = (c^(1-Params.gamma)) / (1-Params.gamma); % utility from consumption
            
            % compute expected value of next period's value function
            EV_next = 0; % initialize expected value of next period's value function
            for k = 1:length(z_grid) % loop over possible next period's income states
                z_next = z_grid(k); % next period's income state
                prob_z_next = z_prob(find(z_grid == z), k); % transition probability to next period's income state
                state_next_idx = find(state_space(:, 1) == a_next & state_space(:, 2) == z_next); % index of next period's state in the state space
                EV_next = EV_next + prob_z_next * V_grid(state_next_idx); % accumulate expected value of next period's value function
            end
            
            expected_values(j) = u + Params.beta * EV_next; % total expected value for current action
        end
        
        [V_grid(i), policy_idx] = max(expected_values); % update value function and policy function for current state
        policy_grid(i) = a_grid(policy_idx); % store optimal next period's wealth as policy function
    end