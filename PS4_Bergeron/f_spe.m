function [v_fxn, policy_fxn, phi_dist, net_assets] = f_spe(r, markov_chain, u_fxn, u_prime_inv, Params, a_grid, a_fine_grid)

% [v_fxn, policy_fxn, phi_dist, net_assets] = 
%   f_spe(r, markov_chain, u_fxn, 
%           u_prime_inv, Params, a_grid, a_fine_grid)

% 

% r = real interestrate
% markov_chan = markov chain with  [Z, PI]
%           Z = state vector
%           PI = transition matrix
% u_fxn: utility function
% u_ptime_inve: inverse of the derivative of the utility function
% a_grid: wealth grid
% a_grid_fine: fine wealth grid for interpolation
% Params: model paramter structures

