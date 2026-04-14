% Econ 527 Spring 2026 Question 1
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 21, 2026
% Last Edited: April 8, 2026


%% Part 1: AR(1) approximation
N = 9; % number of grid points for z
rho = 0.99; % persistence parameter
sigma = 0.2; % standard deviation of the error term
q = 3; % number of standard deviations to cover in the grid

% Tauchen's method
[Z_tauchen, Zprob_tauchen] = tauchen(N, 0, rho, sigma, q);
disp('Tauchen method:');
disp('Grid points:');
disp(Z_tauchen);
disp('Transition probabilities:');
disp(Zprob_tauchen);    

%% Part 2: Simulate histories 

% function to simulate the data 
function X = simluate_markov(Z_probs, n_duration, n_sim)
    
    % calculating the invariant distribution
    invar_dist= mc_invdist(Z_probs);
    disp('Invariant distribution:');
    disp(invar_dist);

    mc = dtmc(Z_probs); % create a discrete-time Markov chain object
    n_states = size(Z_probs, 1);  % use size() directly from Z_probs, not invar_dist


    %figure;
    %graphplot(dtmc_tauchen,'ColorEdges',true);

    % Simulate the Markov chain 1000 times 
    pX = zeros(n_duration + 1, n_sim); % pre-assigning 
    parfor i = 1:n_sim
        
        % Draw initial state from invariant distribution
        x0_idx = find(rand < cumsum(invar_dist), 1);

        % Convert to one-hot vector
        x0 = zeros(1, n_states);
        x0(x0_idx) = 1;

        X(:, i) = simulate(mc, n_duration, 'X0', x0);
    end
end


% function to calculate variables 
function [x_acf, x_var, x_acov1, x_acov2, x_cond_var] = calculate_variables(X, sigma)
    % (a) First-order autocorrelation: cor(xt,xt−1)
    x_acf = autocorr(X(:)); % compute autocorrelation function up to lag 1
    disp('First-order autocorrelation:');
    disp(x_acf(2)); % autocorrelation at lag 1

    % (b) Unconditional variance: var(xt)
    x_var = var(X(:)); % compute variance of the simulated series
    disp('Unconditional variance:');    
    disp(x_var);

    % (c) First-order autocovariance: cov(xt,xt−1)
    x_acov1 = x_acf(2) * x_var; % autocovariance at lag 1 is autocorrelation at lag 1 times variance
    x_acov2 = xcov(X(:), X(:), 1); % compute autocovariance at lag 1 directly
    disp('First-order autocovariance:');
    disp('Calculated from autocorrelation and variance:');
    disp(x_acov1);
    disp('Calculated directly from data:');
    disp(x_acov2);

    % (d) Conditional variance: var(xt|xt−1) = var(εt)
    x_cond_var = sigma^2; % conditional variance is the variance of the error term
    disp('Conditional variance:');
    disp(x_cond_var);   
end

n_duration = 1000; % number of periods to simulate
n_sim = 1000; % number of simulations to run
x_tauchen = simluate_markov(Zprob_tauchen, n_duration, n_sim);

calculate_variables(x_tauchen, sigma);

%% Rouwemnhourt Method

% Rouwenhorst's method
[Z_rouwenhorst, Zprob_rouwenhorst] = rouwenhorst(rho, sigma, N);
disp('Rouwenhorst method:');
disp('Grid points:');
disp(Z_rouwenhorst);
disp('Transition probabilities:');
disp(Zprob_rouwenhorst);  

x_rouwenhorst = simluate_markov(Zprob_rouwenhorst, n_duration, n_sim);
calculate_variables(x_rouwenhorst, sigma);


%% Part 3: Optional 
q = 2.5;
[Z_tauchen_new, Zprob_tauchen_new] = tauchen(N, 0, rho, sigma, q);
disp('Tauchen method with q = 2.5:');
disp('Grid points:');
disp(Z_tauchen_new);
disp('Transition probabilities:');
disp(Zprob_tauchen_new);    

x_tauchen_new = simluate_markov(Zprob_tauchen_new, n_duration, n_sim);
calculate_variables(x_tauchen_new, sigma);

%% Saving results
publish('PS1_Q1.m', 'pdf')