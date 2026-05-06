[V_grid, a_grid_egm] = policy_fxn_egm_calc(Params, z_grid, z_prob, a_prime_grid, V_grid_old, a_grid_egm)   

% [[V_grid, a_grid_egm] = policy_fxn_egm_calc(Params, z_grid, z_prob, a_prime_grid, V_grid_old, a_grid_egm)   
%
% value function iteration loop but for endgenous grid (V^{n+1}) and policy grid (a') for each (a, iz) pair
% takes endogenous grid a' as input and aclaulates
% a(a', iz) for each (a', iz) pair
%
% Params: struct of model parameters (e.g. beta, signa, etc.)
% Z_grid: discretized income grid
% z_prob: transition probabikities for income process
% V_grid_old: current value function (V^n) for each (a, iz) pair
% a_prime_grid: current endogenous wealth grid (a' choices) for each (a, iz)
% a_grid_egm: current exogenous wealth grid (a choices) for each (a', iz)

for ia_prime = 1:Params.n_a_prime
    a_prime = a_prime_grid(ia_prime); % current poly grid
    for iz = 1:Params.n_z
        z= z_grid(iz); % current income 

       % calculate current consumptio
        c = consumption(a_prime, z, a_prime_grid, Params); 
        % consumption for each (a', iz) pair using current endogenous grid a' and income z
        % i.e. next period wealth choice a' and current income z determine current consumption c

        up = utility_prime(c, Params); % marginal utility for each (a', iz) pair

        % calculate expected marginal utility for each (a', iz) pair
        up_exp = 0;
        
        

        % calcualting  
        
