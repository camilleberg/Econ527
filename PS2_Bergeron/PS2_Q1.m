% Econ 527 Spring 2026 Question 2
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 23, 2026
% Last Edited: April 20, 2026

% Neoclassical Growth Model

%% Part 1: Neoclassical Growth Model Setup

% paramters
% for vfi
beta = 0.96; % discount factor
alpha = 0.35; % capital share
delta = 1; % depreciation rate
A = 1; % technology level


%% Anlaytical solution for the value function and policy function
b = alpha / (1 - alpha * beta);
a = (1/(1-beta)) * log((A/(1+beta*b)) * ((A*beta*b)/(1+beta*b))^(beta*b));
disp('Analytical value function: V(k) = a + b*log(k)');
disp(['a = ', num2str(a)]);
disp(['b = ', num2str(b)]);

% creating funcion for the analytical value function
function V = analytical_value_function(k, a, b)
    V = a + b * log(k); % analytical value function
end 

%% Numerical solution using value function iteration

% parameters
% for grid
curve = 2; % curvature parameter for polynomial grid
nodes = 100; % number of grid points
point_a = 0.1; % lower bound of the grid
point_b = 10; % upper bound of the grid

% for iteration 
e_stop = 1e-4; % convergence criterion

% creating grid points for capital using polynomial transformation
function grid_poly = polynomial_grid(point_a, point_b, nodes, curve)
    grid_normalized = linspace(0, 1, nodes); % normalized grid in [0,1]
    grid_poly = point_a + (point_b - point_a) * grid_normalized.^curve; % polynomial transformation
end

% creating capital grid using polynomial transformation
K_grid = polynomial_grid(point_a, point_b, nodes, curve); % capital grid using polynomial transformation
disp('Capital grid points:');
disp(K_grid(1:10));

% cubic spline interpolation for the value function
function V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
    spline_interp = spline(K_grid, V_grid); % create spline interpolation object
    V_interp = ppval(spline_interp, K_query); % evaluate spline interpolation at query points
end 


V_interp = cubic_spline_interpolation(K_grid, V_guess, K_grid); % initial interpolation of the value function at grid points
V_diff = max(abs(V_interp - V_guess)); % initial difference for convergence check
disp('Initial interpolated value function at grid points:');
disp(V_interp(1:10)); % display initial interpolated value function at grid points
disp(['Initial maximum difference for convergence check: ', num2str(V_diff)]); % display initial maximum difference for convergence check


% iteration loop for value function iteration
V_guess = zeros(nodes); % initializing
iteration = 1; % initializing iteration counter
while V_diff > e_stop*(1-beta) % check for convergence
    V_interp = cubic_spline_interpolation(K_grid, V_guess, K_grid); % interpolate value function at grid points
    V_diff = max(abs(V_interp - V_guess)); % calculate maximum difference for convergence check
    disp(['Iteration: ', num2str(iteration), ', Max Difference: ', num2str(V_diff)]); % display iteration and max difference
    V_guess = V_interp; % update value function guess
    iteration = iteration + 1; % increment iteration counter
end 



disp(['Value function iteration converged after ', num2str(iteration), ' iterations.']);
disp('Numerical value function at grid points:');
disp(V_guess(1:10));    
    
% plotting the value function
figure;
fig = figure;
theme(fig, "light");
plot(K_grid, V_guess, 'b-', 'LineWidth', 2); % plot numerical value function
hold on;
K_query = linspace(point_a, point_b, 1000); % query points for plotting
V_analytical = analytical_value_function(K_query, a, b); % analytical value function at query points
plot(K_query, V_analytical, 'r--', 'LineWidth', 2); % plot analytical value function
xlabel('Capital (k)');
ylabel('Value Function (V)');
title('Value Function Iteration vs Analytical Solution');
legend('Numerical Value Function', 'Analytical Value Function');
grid on;
saveas(gcf,'./figs/value_function_comparison.fig')
saveas(gcf,'./figs/value_function_comparison.png')