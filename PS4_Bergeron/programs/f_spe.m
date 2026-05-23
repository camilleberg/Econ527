function [v_fxn, policy_fxn, phi_dist, net_assets] = f_spe(z_grid, z_prob, u_fxn, u_prime_inv, Params, a_next_grid, a_fine_grid)

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

%% ── Setup ────────────────────────────────────────────────────────────────

R = 1 + Params.r;

assert(isscalar(Params.r), 'Params.r must be a scalar; got size %dx%d', size(Params.r,1), size(Params.r,2));
assert(isscalar(R),        'R must be a scalar');

% One-argument closures so subfunctions receive pre-baked Params
u_fxn_h       = @(c)  u_fxn(c,       Params);
u_prime_inv_h = @(mu) u_prime_inv(mu, Params);

a_1 = a_next_grid(1);   % borrowing limit = first grid point

%% ── Initialisation ───────────────────────────────────────────────────────

% Initial consumption: spend everything above the borrowing limit
c_init = zeros(Params.n_a, Params.n_z);
for iz = 1:Params.n_z
    c_init(:, iz) = max(R * a_next_grid + exp(z_grid(iz)) - a_1, 1e-6);
end

% Initial marginal utility: u'(c) = c^{-gamma}
mu_grid = c_init .^ (-Params.gamma);

%% ── EGM loop — iterate on mu only, no V update inside ───────────────────
%
% Key design decision: updating V inside the loop via spline interpolation
% on the endogenous grid causes explosive extrapolation errors because the
% endogenous grid points do not cover the full exogenous grid range each
% iteration. Instead we iterate only on mu (equivalently c), which is
% well-behaved, and recover V once after convergence.

error_val = Inf;
iteration = 0;
tic;

while error_val > Params.e_stop && iteration < Params.max_iter

    mu_grid_old = mu_grid;

    %── Step 1: EGM backward step ─────────────────────────────────────────
    % Given mu(a',z), invert Euler eq to get c(a',z) and endogenous a(a',z)
    [a_current_grid, c_current_grid] = compute_EGM( ...
        Params, z_grid, z_prob, a_next_grid, mu_grid_old, u_prime_inv_h);

    %── Step 2: classify constrained / unconstrained ──────────────────────
    % Unconstrained: endogenous a > a_min (FOC holds with equality)
    % Constrained:   endogenous a <= a_min (household would borrow more
    %                than allowed; they are forced to save exactly a_1)
    constrained_mask   = a_current_grid <= Params.a_min;
    unconstrained_mask = ~constrained_mask;

    %── Step 3: build c_new on the exogenous grid ─────────────────────────
    c_new = zeros(Params.n_a, Params.n_z);

    for iz = 1:Params.n_z
        unc = unconstrained_mask(:, iz);
        con = constrained_mask(:, iz);

        %── Unconstrained: interpolate c from endogenous → exogenous grid ─
        % Sort and deduplicate endogenous grid before splining
        a_endog = real(a_current_grid(unc, iz));
        c_endog = real(c_current_grid(unc, iz));
        [a_s, si] = sort(a_endog);
        c_s       = c_endog(si);
        [a_u, ui] = unique(a_s);
        c_u       = c_s(ui);

        if length(a_u) >= 2
            query   = a_next_grid(unc);
            unc_idx = find(unc);

            % ONLY interpolate — never extrapolate (lecture note requirement)
            % Points outside the endogenous range are handled as constrained
            in_range  = query >= a_u(1) & query <= a_u(end);
            out_range = ~in_range;

            if any(in_range)
                c_new(unc_idx(in_range), iz) = max( ...
                    cubic_spline_interpolation(a_u, c_u, query(in_range)), 1e-10);
            end
            % Out-of-range unconstrained points: use constrained formula
            % (they sit below the lowest endogenous gridpoint, so the
            %  borrowing constraint is effectively binding there too)
            if any(out_range)
                a_out = a_next_grid(unc_idx(out_range));
                c_new(unc_idx(out_range), iz) = max( ...
                    R * a_out + exp(z_grid(iz)) - a_1, 1e-10);
            end
        else
            % Too few unique endogenous points: fall back to constrained formula
            unc_idx = find(unc);
            c_new(unc_idx, iz) = max( ...
                R * a_next_grid(unc_idx) + exp(z_grid(iz)) - a_1, 1e-10);
        end

        %── Constrained: consume cash-on-hand minus minimum saving ────────
        c_new(con, iz) = max( ...
            R * a_next_grid(con) + exp(z_grid(iz)) - a_1, 1e-10);
    end

    %── Step 4: update mu from new consumption ────────────────────────────
    mu_grid = c_new .^ (-Params.gamma);

    %── Step 5: convergence on relative change in mu ──────────────────────
    % Relative criterion avoids scale problems (mu near constraint is huge)
    error_val = max(abs(mu_grid(:) - mu_grid_old(:)) ./ (abs(mu_grid_old(:)) + 1e-10));
    iteration = iteration + 1;

    if mod(iteration, 500) == 0
        fprintf('  EGM iter %d, rel error = %.2e\n', iteration, error_val);
    end
end
time_egm = toc;

fprintf('EGM converged after %d iterations in %.2f seconds (error = %.2e)\n', ...
        iteration, time_egm, error_val);

%% ── Recover policy grid from converged consumption ───────────────────────
% Budget constraint: R*a + exp(z) = c + a'  =>  a' = R*a + exp(z) - c

policy_grid = zeros(Params.n_a, Params.n_z);
for iz = 1:Params.n_z
    policy_grid(:, iz) = R * a_next_grid + exp(z_grid(iz)) - c_new(:, iz);
    policy_grid(:, iz) = max(policy_grid(:, iz), a_1);   % enforce constraint
end

%% ── Compute V once via Howard iteration with fixed policy ────────────────
% Now that the policy is converged, V satisfies the fixed-policy Bellman:
%   V(a,z) = u(R*a + exp(z) - g(a,z)) + beta * E[V(g(a,z), z') | z]
% This is a linear equation in V — Howard iteration converges fast.

% Initial guess: perpetual consumption discounted forever
V_grid = zeros(Params.n_a, Params.n_z);
for iz = 1:Params.n_z
    c_pol = max(R * a_next_grid + exp(z_grid(iz)) - policy_grid(:, iz), 1e-10);
    V_grid(:, iz) = u_fxn_h(c_pol) / (1 - Params.beta);
end

% Howard iteration (linear — no maximisation needed)
for howard_iter = 1:2000
    V_new = zeros(Params.n_a, Params.n_z);
    for iz = 1:Params.n_z
        % Expected continuation value: E[V(g(a,z), z') | z]
        EV = zeros(Params.n_a, 1);
        for iz2 = 1:Params.n_z
            % Interpolate V at next-period assets chosen by policy
            a_pol_iz = max(a_next_grid(1), min(a_next_grid(end), policy_grid(:, iz)));
            EV = EV + z_prob(iz, iz2) * ...
                cubic_spline_interpolation(a_next_grid, V_grid(:, iz2), a_pol_iz);
        end
        c_pol = max(R * a_next_grid + exp(z_grid(iz)) - policy_grid(:, iz), 1e-10);
        V_new(:, iz) = u_fxn_h(c_pol) + Params.beta * EV;
    end

    v_err   = max(abs(V_new(:) - V_grid(:)));
    V_grid  = V_new;
    if v_err < 1e-8
        fprintf('Howard iteration converged after %d steps (V error = %.2e)\n', ...
                howard_iter, v_err);
        break
    end
end

%% ── Interpolate value and policy onto fine grid ──────────────────────────

n_a_fine   = length(a_fine_grid);
v_fxn      = zeros(n_a_fine, Params.n_z);
policy_fxn = zeros(n_a_fine, Params.n_z);

for iz = 1:Params.n_z
    v_fxn(:, iz)      = cubic_spline_interpolation( ...
        a_next_grid, V_grid(:, iz),      a_fine_grid);
    policy_fxn(:, iz) = cubic_spline_interpolation( ...
        a_next_grid, policy_grid(:, iz), a_fine_grid);
    % Enforce borrowing constraint on interpolated policy
    policy_fxn(:, iz) = max(policy_fxn(:, iz), a_fine_grid(1));
end

%% ── Build state-transition matrix on fine grid ───────────────────────────
% policy_P(s', s) = prob of transitioning from state s to state s'
% States ordered as: s = (iz-1)*n_a_fine + ia

numstates = n_a_fine * Params.n_z;
policy_P  = zeros(numstates);

for iz = 1:Params.n_z
    for ia = 1:n_a_fine
        s_from  = (iz - 1) * n_a_fine + ia;

        a_prime = policy_fxn(ia, iz);
        a_prime = max(a_fine_grid(1), min(a_fine_grid(end), a_prime));

        % Linear interpolation bracket
        ik  = min(find(a_fine_grid <= a_prime, 1, 'last'), n_a_fine - 1);
        w_lo = (a_fine_grid(ik+1) - a_prime)        / (a_fine_grid(ik+1) - a_fine_grid(ik));
        w_hi = (a_prime           - a_fine_grid(ik)) / (a_fine_grid(ik+1) - a_fine_grid(ik));

        for iz_next = 1:Params.n_z
            pi_z  = z_prob(iz, iz_next);
            s_lo  = (iz_next - 1) * n_a_fine + ik;
            s_hi  = (iz_next - 1) * n_a_fine + (ik + 1);
            policy_P(s_lo, s_from) = policy_P(s_lo, s_from) + pi_z * w_lo;
            policy_P(s_hi, s_from) = policy_P(s_hi, s_from) + pi_z * w_hi;
        end
    end
end

%% ── Stationary distribution via eigenvalue method ────────────────────────
% mc_invdist issues a warning if the unit eigenvalue is not unique

stationary_eigen = mc_invdist(policy_P');   % (numstates x 1), sums to 1

% Reshape to (n_a_fine x n_z) joint distribution
phi_dist = reshape(stationary_eigen, n_a_fine, Params.n_z);

%% ── Net asset demand ─────────────────────────────────────────────────────
% E[a] = sum_i sum_j a_i * phi(a_i, z_j)

lambda_marginal_a = sum(phi_dist, 2);                          % (n_a_fine x 1)
net_assets        = sum(a_fine_grid(:) .* lambda_marginal_a);  % scalar

end


%% ══════════════════════════════════════════════════════════════════════════
%  Local subfunctions
%% ══════════════════════════════════════════════════════════════════════════

% ── cubic_spline_interpolation ────────────────────────────────────────────
function V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
% V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
%
% Cubic spline interpolation. Non-finite values in V_grid are replaced with
% a large negative number before splining to avoid NaN propagation.
%
% NOTE: this function CAN extrapolate if K_query falls outside K_grid.
%       The caller is responsible for restricting K_query to the interior
%       when extrapolation would be invalid (e.g. in the EGM c update).

    V_grid_clean = V_grid;
    V_grid_clean(~isfinite(V_grid_clean)) = -1e10;

    spline_interp = spline(K_grid, V_grid_clean);
    V_interp      = ppval(spline_interp, K_query);
end


% ── compute_EGM ──────────────────────────────────────────────────────────
function [a_current_grid, c_current_grid] = compute_EGM( ...
        Params, z_grid, z_prob, a_next_grid, mu_grid, u_prime_inv)
% [a_current_grid, c_current_grid] = compute_EGM(...)
%
% EGM backward step. Given mu(a',z) = u'(c(a',z)) on the exogenous grid,
% inverts the Euler equation to recover c and the endogenous wealth grid a:
%
%   Euler:  u'(c) = beta*(1+r) * E[mu(a',z') | z]
%   =>  c = (u')^{-1}( beta*(1+r) * E[mu] )
%   =>  a = (c + a' - exp(z)) / (1+r)          [from budget constraint]
%
% Inputs:
%   mu_grid     - (n_a x n_z) current marginal utility on exogenous grid
%   u_prime_inv - one-argument handle: (u')^{-1}(x) = x^{-1/gamma}

    beta = Params.beta;
    R    = 1 + Params.r;

    % E[mu(a'_i, z') | z_j] = sum_{z'} pi(z_j, z') * mu(a'_i, z')
    % mu_grid: (n_a x n_z),  z_prob: (n_z x n_z)
    % EMU = mu_grid * z_prob'  gives (n_a x n_z)  [correct orientation]
    EMU = mu_grid * z_prob';

    a_current_grid = zeros(Params.n_a, Params.n_z);
    c_current_grid = zeros(Params.n_a, Params.n_z);

    for ia = 1:Params.n_a
        a_next = a_next_grid(ia);
        for iz = 1:Params.n_z
            z_curr = z_grid(iz);

            % Invert Euler: c = (u')^{-1}(beta*R * E[mu])
            c_curr = max(u_prime_inv(beta * R * EMU(ia, iz)), 1e-10);

            % Recover current wealth from budget constraint:
            %   R*a + exp(z) = c + a'  =>  a = (c + a' - exp(z)) / R
            a_curr = (c_curr + a_next - exp(z_curr)) / R;

            a_current_grid(ia, iz) = a_curr;
            c_current_grid(ia, iz) = c_curr;
        end
    end
end


% ── mc_invdist ───────────────────────────────────────────────────────────
function P = mc_invdist(PI)
% P = mc_invdist(PI)
%
% Computes the invariant distribution of a Markov chain with row-stochastic
% transition matrix PI by finding the unit left eigenvector of PI.

    [V, D] = eig(PI');

    ii = find(abs(diag(D) - 1) < 1e-8, 1);   % index of unit eigenvalue

    P = V(:, ii) / sum(V(:, ii));             % normalise to sum to 1

    assert(max(abs(P' - P' * PI)) < 1e-12, ...
           'mc_invdist: stationary distribution verification failed');

    if sum(abs(diag(D) - 1) < 1e-8) > 1
        warning('mc_invdist: unit eigenvalue not unique — distribution may not be unique');
    end
end