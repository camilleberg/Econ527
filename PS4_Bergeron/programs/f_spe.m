function [V_grid, policy_fxn, phi_dist, net_assets] = f_spe(r, z_grid, z_prob, u_fxn, u_prime_inv, Params, a_grid, a_fine_grid)

% [v_fxn, policy_fxn, phi_dist, net_assets] =
%       f_spe(z_grid, z_prob, u_fxn, u_prime_inv, Params, a_next_grid, a_fine_grid)
%
% Computes the stationary partial equilibrium of a Huggett / Bewley model
% using the Endogenous Grid Method (EGM) and returns:
%
%   v_fxn      - (n_a_fine x n_z)  value function on fine grid
%   policy_fxn - (n_a_fine x n_z)  savings policy function on fine grid
%   phi_dist   - (n_a_fine x n_z)  stationary joint distribution phi(a,z)
%   net_assets - scalar, sum_i sum_j a_i * phi(a_i, z_j)
%
% Algorithm:
%   1. EGM loop iterates on marginal utility mu ONLY (not V). This avoids
%      the numerical blowup caused by spline extrapolation in V updates.
%   2. Convergence is checked on relative change in mu (cleaner scale).
%   3. After EGM converges, policy_grid is recovered from the budget
%      constraint. V is then computed ONCE via Howard iteration with the
%      fixed converged policy — fast and numerically stable.
%
% Inputs:
%   z_grid       - (n_z x 1) log-income grid
%   z_prob       - (n_z x n_z) transition matrix (rows sum to 1)
%   u_fxn        - utility handle: u_fxn(c, Params)
%   u_prime_inv  - inverse marginal utility handle: u_prime_inv(mu, Params)
%   Params       - struct with fields:
%                    .r, .beta, .gamma
%                    .a_min, .a_max, .n_a, .n_z, .curve
%                    .max_iter, .e_stop
%   a_next_grid  - (n_a x 1) exogenous savings grid
%   a_fine_grid  - (n_a_fine x 1) fine grid for output

%% initialize
% One-argument closures so subfunctions receive pre-baked Params


V_grid = NaN(Params.n_z,Params.n_a); % initilaize value function
dV0_grid = NaN(Params.n_z,Params.n_a); % intialize derivative at value function

policy_fxn = NaN(Params.n_z,Params.n_af);

pp = cell(Params.n_z,1); % for spline interpolanrs
pp_d1 = cell(Params.n_z,1); % for 


%% 1. solve bellman

c_init = max(r,0.001) * (a_grid-a_grid(1)+0.01) + z_grid; % initial consumption 
V_init = u_fxn(c_init) / (1-Params.beta);    % Initial guess for V


% iteration loop 
for iter = 1:Params.max_iter
    % looping through income states and interpolating
    for iz = 1:Params.n_z
        pp{iz} = spline(a_grid, V_init(iz,:));
        pp_d1{iz} = fnder(pp{iz}, 1);
        dV0_grid(iz,:) = ppval(pp_d1{iz}, a_grid);
    end

    % calculating the continuation value 
    EV = z_prob * V_init; % expected income state 
    EDV = z_prob * dV0_grid; % expected marginal

    % caluclate endogeneous grid points a 
    % clamp EDV to be positive before passing to u_prime_inv to avoid complex values
    a_grid_end = (a_grid - z_grid + u_prime_inv(Params.beta*max(EDV, 1e-10))) / (1+r);

    % update the value function for endogenous grids 
    V_grid_end = u_fxn((1+r)*a_grid_end + z_grid - a_grid) + Params.beta*EV;

    % calculating contrainsed and unconstrained 
    for iz = 1:Params.n_z
        % updating unconsrained grid / interpolating regularly 
      V_grid(iz,:) = spline(a_grid_end(iz,:), V_grid_end(iz,:), a_grid);    
      
      % calculating constrained mask to find exogeneous grid
      constrained_mask = a_grid < a_grid_end(iz,1);   
      
      % updating value function at constrained point to be the budget constrain
      V_grid(iz,constrained_mask) = u_fxn((1+r)*a_grid(constrained_mask) + z_grid(iz) - ...
        a_grid(1)) + Params.beta*EV(iz,1);
   end

   % Check for convergence — use (1-beta)-scaled tolerance to match level of V
   tol = (1 - Params.beta) * 1e-4;
   error = max(max(abs(V_grid - V_init)));  % disp(dif)
   
   if error < tol
      break
   end
   
   % reassign 
   V_init = V_grid;
end

assert(iter < Params.max_iter)

%% 2. Compute policy function on fine wealth grid

% looping over income states
for iz = 1:Params.n_z    
    % interpolating on finer welath grid
   policy_fxn(iz,:) = spline(a_grid_end(iz,:), a_grid, a_fine_grid);
   constrained_mask = a_fine_grid < a_grid_end(iz,1);
    
   % changinf the policy grid to the borrowing limit when hh are constraines
   policy_fxn(iz,constrained_mask) = a_grid(1);

end

assert(min(min(policy_fxn)) >= a_grid(1))

% impose that a' <= a_max (i.e. stays inside welath grid)
policy_fxn = min(policy_fxn, a_grid(Params.n_a));    

%% 3. compute transition matrix

transition_matrix = zeros(Params.n_af * Params.n_z);

% looping through all wealth states 
for ia = 1:Params.n_af
    % looping through income states 
    for iz = 1:Params.n_z
    
        j = min(sum(a_fine_grid <= policy_fxn(iz,ia)), Params.n_af-1); 
        wgt = (policy_fxn(iz,ia) - a_fine_grid(j)) / (a_fine_grid(j+1) - a_fine_grid(j));  
    
        assert(0<=wgt && wgt<=1)
        s = (iz-1)*Params.n_af + ia;
        
        % looping through income again
        for iz_new = 1:Params.n_z
            s1_lo = (iz_new-1)*Params.n_af + j;
            transition_matrix(s,s1_lo) = (1-wgt) * z_prob(iz,iz_new); 
            transition_matrix(s,s1_lo+1) = wgt * z_prob(iz,iz_new); 
            
        end
        
    end
end

assert(max(abs(sum(transition_matrix,2) - 1)) < 1E-12)

%% 4. calculate invariant distribution 
phi_vec = mc_invdist(transition_matrix); % using eigen values
phi_dist = reshape(phi_vec,Params.n_af,Params.n_z)';  % reshaping to use            

% calculating net asset demand
net_assets = sum(sum(a_fine_grid.*phi_dist));          

end

%% local functions

