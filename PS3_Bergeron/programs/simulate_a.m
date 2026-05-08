function a_sim = simulate_a(policy_fxn, a_grid, z_sim, Params)
% a_sim = simulate_a(policy_fxn, a_grid, z_sim, Params)
%
% policy_fxn:  policy function indices (n_assets x n_states)
% a_grid:      asset grid
% z_sim:       simulated income states (n_duration+1 x n_sim)


    n_duration  = Params.n_duration;
    n_sim       = Params.n_sim;

    % convert policy function values to indices on a_grid
    policy_idx = zeros(size(policy_fxn));
    
    % tried to use discrete index earlier but that didn;t work
    % policy fxn is continuous, so need to interpoate 
    for iz = 1:size(policy_fxn, 2)
        policy_idx(:, iz) = interp1(a_grid, 1:length(a_grid), ...
                                    policy_fxn(:, iz), 'nearest', 'extrap');
    end
    policy_idx = round(policy_idx);   % ensure integer indices

    % find index of borrowing limit
    a0_idx = 1;   % a_min is first point on grid by construction

    a_idx_sim       = zeros(n_duration + 1, n_sim);
    a_idx_sim(1, :) = a0_idx;

    for t = 1:n_duration
        a_idx_sim(t+1, :) = policy_idx(sub2ind(size(policy_idx), ...
                                a_idx_sim(t,:), z_sim(t,:)));
    end

    a_sim = a_grid(a_idx_sim);
end