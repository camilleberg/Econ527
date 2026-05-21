function net_assets = f_nad(z_grid, z_prob, u_fxn, u_prime_inv, Params, a_grid, a_fine_grid)
% net_assets = f_nad(r, markov_chain, u_fxn, u_prime_inv, Params, a_grid, a_fine_grid)
%
% Wrapper around f_spe that returns only net asset demand.
% Useful for root-finding routines (e.g. fzero) that require scalar output.
%
% Inputs:  identical to f_spe — see f_spe.m for full documentation
% Output:
%   net_assets - scalar, sum_i sum_j a_i * phi(a_i, z_j)
%                = 0 at a general equilibrium interest rate

    [~, ~, ~, net_assets] = f_spe(z_grid, z_prob, @u_fxn, @u_prime_inv, Params, a_grid, a_fine_grid);
end