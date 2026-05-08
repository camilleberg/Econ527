function [a_current_grid, c_current_grid] = compute_EGM(Params, z_grid, z_prob, a_next_grid, mu_grid)
 
% [a_current_grid, c_current_grid] = compute_EGM(Params, z_grid, z_prob, a_next_grid, mu_grid)
%
% EGM backward step. Given the exogenous grid a' and the current marginal
% utility grid mu(a',z) = u'(c(a',z)) = c(a',z)^{-gamma}, recovers the
% endogenous current-period wealth grid a(a',z) by inverting the Euler eq:
%
%   Euler:  u'(c) = beta*(1+r)*E[u'(c(a',z'))|z]
%   =>  c(a',z) = [beta*(1+r)*E[mu(a',z')|z]]^{-1/gamma}
%   =>  a(a',z) = (a' - exp(z) + c) / (1+r)
%
% Inputs:
%   Params      - parameter struct
%   z_grid      - (n_z x 1) discretised income grid
%   z_prob      - (n_z x n_z) transition matrix
%   a_next_grid - (n_a x 1) exogenous savings grid (fixed throughout)
%   mu_grid     - (n_a x n_z) marginal utility u'(c(a',z)) from last iter
%
% Outputs:
%   a_current_grid - (n_a x n_z) endogenous current-wealth grid
%   c_current_grid - (n_a x n_z) consumption implied by Euler equation
 
    beta = Params.beta;
    R    = 1 + Params.r;
 
    % BUG FIX 1: expected marginal utility via transition matrix
    % EMU(ia,iz) = sum_{iz'} pi(iz,iz') * mu(ia,iz')
    % mu_grid is (n_a x n_z), z_prob is (n_z x n_z)
    % We want row iz of z_prob dotted with columns of mu_grid
    % => EMU = mu_grid * z_prob'   (n_a x n_z)
    EMU = mu_grid * z_prob';
 
    a_current_grid = zeros(Params.n_a, Params.n_z);
    c_current_grid = zeros(Params.n_a, Params.n_z);
 
    for ia = 1:Params.n_a
        a_next = a_next_grid(ia);
        for iz = 1:Params.n_z
            z_curr = z_grid(iz);
 
            % BUG FIX 2: Euler equation requires beta*R factor (not just beta)
            % u'(c) = beta*(1+r)*E[u'(c')|z]  =>  c = (beta*R*EMU)^{-1/gamma}
            rhs    = max(beta * R * EMU(ia, iz), 1e-10);  % clamp to keep real
            c_curr = rhs ^ (-1 / Params.gamma);
 
            % BUG FIX 3: was using gradient(V_grid_old) which is wrong —
            % EGM tracks marginal utility directly, no gradient of V needed.
            % Budget: c = R*a + exp(z) - a'  =>  a = (a' - exp(z) + c) / R
            % (matches problem set formula: a = (a' - exp(z) + u'^{-1}[...])/R)
            a_curr = (a_next - exp(z_curr) + c_curr) / R;
 
            a_current_grid(ia, iz) = a_curr;
            c_current_grid(ia, iz) = c_curr;
        end
    end
end
 