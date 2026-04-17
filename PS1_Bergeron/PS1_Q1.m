% Econ 527 Spring 2026 Question 1
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 21, 2026
% Last Edited: April 8, 2026

%% Part 1: AR(1) approximation
N = 9;      % number of grid points for z
rho = 0.99; % persistence parameter
sigma = 0.2; % standard deviation of the error term
q = 3;      % number of standard deviations to cover in the grid

% Tauchen's method
[Z_tauchen, Zprob_tauchen] = tauchen(N, 0, rho, sigma, q);
disp('Tauchen method:');
disp('Grid points:'); disp(Z_tauchen);
disp('Transition probabilities:'); disp(Zprob_tauchen);

%% Part 2: Simulate histories

function X = simulate_markov(Z_probs, n_duration, n_sim)

    % Invariant distribution
    invar_dist = mc_invdist(Z_probs);
    disp('Invariant distribution:'); disp(invar_dist);

    mc = dtmc(Z_probs);
    n_states = size(Z_probs, 1);

    % Plotting — create figure first, then plot into it
    %fig = figure;
    %theme(fig, "light");
    %graphplot(mc, 'ColorEdges', true);
    %varName = inputname(1);
    %saveas(fig, varName + ".png");

    % Simulate the Markov chain
    X = zeros(n_duration + 1, n_sim); % pre-allocate X (not pX)
    parfor i = 1:n_sim
        x0_idx = find(rand < cumsum(invar_dist), 1);
        x0 = zeros(1, n_states);
        x0(x0_idx) = 1;
        X(:, i) = simulate(mc, n_duration, 'X0', x0);
    end
end

function [x_acf, x_var, x_acov1, x_cond_var] = calculate_variables(X, sigma)
    % (a) First-order autocorrelation
    x_autocor = autocorr(X(:));
    x_acf = x_autocor(2);
    disp('First-order autocorrelation:'); disp(x_acf);

    % (b) Unconditional variance
    x_var = var(X(:));
    disp('Unconditional variance:'); disp(x_var);


    % (c) First-order autocovariance
    x_acov1 = x_acf * x_var;
    disp('First-order autocovariance:');
    disp('From autocorrelation * variance:'); disp(x_acov1);

    % (d) Conditional variance
    x_cond_var = sigma^2;
    disp('Conditional variance:'); disp(x_cond_var);
end

n_duration = 1000;
n_sim = 1000;

x_tauchen = simulate_markov(Zprob_tauchen, n_duration, n_sim);
vars_tauchen = calculate_variables(x_tauchen, sigma);

%% Rouwenhorst Method
[Z_rouwenhorst, Zprob_rouwenhorst] = rouwenhorst(rho, sigma, N);
disp('Rouwenhorst method:');
disp('Grid points:'); disp(Z_rouwenhorst);
disp('Transition probabilities:'); disp(Zprob_rouwenhorst);

x_rouwenhorst = simulate_markov(Zprob_rouwenhorst, n_duration, n_sim);
vars_rouwenhorst = calculate_variables(x_rouwenhorst, sigma);

%% Part 3: Optional
q = 2.5;
[Z_tauchen_new, Zprob_tauchen_new] = tauchen(N, 0, rho, sigma, q);
disp('Tauchen method with q = 2.5:');
disp('Grid points:'); disp(Z_tauchen_new);
disp('Transition probabilities:'); disp(Zprob_tauchen_new);

x_tauchen_new = simulate_markov(Zprob_tauchen_new, n_duration, n_sim);
vars_tauchen_new = calculate_variables(x_tauchen_new, sigma);

%% Theoretical values
actual_acf      = rho;
actual_var      = sigma^2 / (1 - rho^2);
actual_acov1    = rho * actual_var;
actual_cond_var = sigma^2;

vars_actual = [actual_acf, actual_var, actual_acov1, actual_acov2, actual_cond_var];
disp('Theoretical values:'); disp(vars_actual);

disp('Tauchen as ratio to theoretical:');
disp(vars_tauchen ./ vars_actual);

disp('Rouwenhorst as ratio to theoretical:');
disp(vars_rouwenhorst ./ vars_actual);  