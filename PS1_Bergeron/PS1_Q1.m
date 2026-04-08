% Econ 527 Spring 2026
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 21, 2026
% Last Edited: April 8, 2026

%% Part 1: AR(1) approximation
n_z = 9; % number of grid points for z
rho = 0.99; % persistence parameter
sigma = 0.2; % standard deviation of the error term
q = 3; % number of standard deviations to cover in the grid
z_grid = linspace(-q*sigma, q*sigma, n_z); % grid for z w q st devs
disp('Grid for z:');
disp(z_grid);

%% Part 2: Simulate histories 
P = dtmc(z_grid); 
disp(P)
T = 1000; % number of time periods
z_sim = zeros(T, 1); % pre-allocate for z

X = simulate(mc,numSteps); % simulate Markov chain to get indices of z_grid
z_sim = simulate_ar1(rho, sigma, T); % simulate AR(1) process




