% Econ 527 Spring 2026
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

% calculating the invariant distribution
invar_dist= mc_invdist(Zprob_tauchen);
disp('Invariant distribution:');
disp(invar_dist);

n_duration = 1000; % number of periods to simulate
n_sim = 1000; % number of simulations
dtmc_tauchen = dtmc(Zprob_tauchen); % create a discrete-time Markov chain object

%figure;
%graphplot(dtmc_tauchen,'ColorEdges',true);

% Simulate the Markov chain 1000 times 
parfor i = 1:n_sim
    xo = find(rand < cumsum(invar_dist), 1); % initial state drawn from invariant distribution
    X(:, i) = simulate(dtmc_tauchen, n_duration);
end

% figure;
% simplot(dtmc_tauchen,X)

% checkin mean value 
disp('Mean value of simulated series:');
disp(mean(X(:)));


% (a) First-order autocorrelation: cor(xt,xt−1)
x_acf = autocorr(X(:)); % compute autocorrelation function up to lag 1
disp('First-order autocorrelation:');
disp(x_acf(2)); % autocorrelation at lag 1

% (b) Unconditional variance: var(xt)
% (c) First-order autocovariance: cov(xt,xt−1)
% (d) Conditional variance: var(xt|xt−1) = var(εt)








