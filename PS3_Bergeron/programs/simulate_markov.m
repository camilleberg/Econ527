function markov_sim = simulate_markov(P, n_duration, n_sim)

% markov_sim = simulate_markov(P, n_duration, n_sim)
% 
% calculates markov chain 
% 
% P: transition state matrix
% n_duration: length of simulation
% n_sim: bnumber of simulations
% invar_dist: statioanry distribution to draw from


    % Invariant distribution
    invar_dist = mc_invdist(P);
    disp('Invariant distribution:'); disp(invar_dist);

    mc = dtmc(P);
    n_states = size(P, 1);

    % Simulate the Markov chain
    markov_sim = zeros(n_duration + 1, n_sim); % pre-allocate X (not pX)
    parfor i = 1:n_sim
        % draws from distribution
        x0_idx = find(rand < cumsum(invar_dist), 1);

        % intial state
        x0 = zeros(1, n_states);

        % setting equal to 1, but now changed to borrowing limit
        x0(x0_idx) = 1;
        markov_sim(:, i) = simulate(mc, n_duration, 'X0', x0);
    end
end