% Econ 527 Spring 2026 
% HW 2 Question 1
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 30, 2026
% Last Edited: April 20, 2026


% Neoclassical Growth Model

%% Part 1: Neoclassical Growth Model Setup

% paramters
% for vfi
Params.beta = 0.96; % discount factor
Params.alpha = 0.35; % capital share
Params.delta = 1; % depreciation rate
Params.A = 1; % technology level


%% Anlaytical solution for the value function and policy function
b = Params.alpha / (1 - Params.alpha * Params.beta);
a = (1/(1-Params.beta)) * log((Params.A/(1+Params.beta*b)) * ((Params.A*Params.beta*b)/(1+Params.beta*b))^(Params.beta*b));
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

% creating capital grid using polynomial transformation
K_grid = polynomial_grid(point_a, point_b, nodes, curve); % capital grid using polynomial transformation
disp('Capital grid points:');
disp(K_grid(1:10));

% cubic spline interpolation for the value function
function V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
    spline_interp = spline(K_grid, V_grid); % create spline interpolation object
    V_interp = ppval(spline_interp, K_query); % evaluate spline interpolation at query points
end 

% value function iteration
V_grid = zeros(size(K_grid)); % initialize value function grid
policy_grid = zeros(size(K_grid)); % initialize policy function grid
iteration = 0; % iteration counter
error = Inf; % initialize error

% running through loop 
while error > e_stop * (1-Params.beta) % continue until convergence
    V_grid_old = V_grid; % store old value function grid
    for i = 1:length(K_grid)
        k = K_grid(i); % current capital stock
        % compute the value of consuming all output and investing the rest
        consumption = Params.A * k^Params.alpha - K_grid; % consumption for each possible next period capital stock
        utility = log(consumption); % utility from consumption
        V_next = cubic_spline_interpolation(K_grid, V_grid_old, K_grid); % interpolate value function for next period capital stocks
        total_value = utility + Params.beta * V_next; % total value for each possible next period capital stock
        [V_grid(i), policy_index] = max(total_value); % update value function and policy function
        policy_grid(i) = K_grid(policy_index); % store optimal next period capital stock in policy grid
    end
    error = max(abs(V_grid - V_grid_old)); % compute maximum error between iterations
    iteration = iteration + 1; % increment iteration counter
end

disp(['Value function iteration converged in ', num2str(iteration), ' iterations.']);
disp('Optimal policy function (next period capital stock):');
disp(policy_grid(1:10));

%% Plotting the results
% Plot value function
figure;
fig = figure;
theme(fig, "light");
plot(K_grid, V_grid, 'b-', 'LineWidth', 2);
hold on;
plot(K_grid, analytical_value_function(K_grid, a, b), 'r--', 'LineWidth', 2);
xlabel('Capital Stock (k)');
ylabel('Value Function (V)');
title('Value Function Iteration vs Analytical Solution');
legend('Numerical VFI', 'Analytical Solution');
grid on;    
saveas(gcf, 'figs/Value_Function_Comparison.png');