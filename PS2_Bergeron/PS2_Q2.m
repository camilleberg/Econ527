% Econ 527 Spring 2026 
% HW 2 Question 2
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 30, 2026
% Last Edited: April 20, 2026


% Income Fluctuations

%% Setting up Problem 

tic; 

% creating parameters
Params.rho = 0.9; % persistence of income process
Params.sigma = 0.2; % standard deviation of income shocks
Params.beta = 0.96; % discount factor
Params.gamma = 1.5; % risk aversion parameter
Params.r = 0.02; % interest rate
Params.a_min = 0; % minimum income level
Params.a_max = 50; % maximum income level
Params.n_a = 100; % number of grid points for wealth grid
Params.curve = 2; % curvature parameter for polynomial grid
Params.n_z = 5; % number of points in discrete approximation for income process
Params.max_iter = 1000; % maximum number of iterations for value function iteration 
Params.e_stop = 1e-4; % convergence criterion for value function iteration


% creating polynomial wealth grid 
a_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve); % income grid using polynomial transformation
disp('Wealth grid points:');
disp(a_grid(1:10));

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method
disp('Transition probabilities:'); disp(z_prob);


%% Value Function Iteration

% initialize value function and policy function
V_grid = zeros(Params.n_a, Params.n_z); % initialize value function grid
policy_grid = zeros(Params.n_a, Params.n_z); % initialize policy function grid
iteration = 0; % initialize iteration counter
error = Inf; % initialize error for convergence check

while error > Params.e_stop * (1-Params.beta) && iteration < Params.max_iter % continue until convergence or max iterations

    V_grid_old = V_grid; % store old value function grid for convergence check

    for ia = 1:Params.n_a
        a = a_grid(ia);
        for iz = 1:Params.n_z
            z = z_grid(iz);

             % consumption for each possible next-period wealth choice
             c_vec = consumption(a, z, a_grid, Params); % (Na x 1)

             % utility 
             u = utility(c_vec, Params); % calcualting current utility

             % continuation
             EV = zeros(Params.n_a, 1); % NA x 1

             % looping over all possible next income values 
            for iz_next = 1:Params.n_z
                % interpolate V(:, iz_next) at a_grid query points
                V_next = cubic_spline_interpolation(a_grid, V_grid_old(:, iz_next), a_grid);
                EV = EV + z_prob(iz, iz_next) * V_next; % (Na x 1)
            end

             total_val = u + Params.beta * EV; % NA x 1 
             total_val = total_val(:);              % force column vector to be safe

             % optimizing
             %[V_grid(ia, iz), policy_grid(ia, iz)] = max(total_val);

             [V_grid(ia, iz), best_idx] = max(total_val);
             if isempty(best_idx) || best_idx < 1 || best_idx > Params.n_a
                best_idx = 1;  % fallback to minimum savings
             end
             policy_grid(ia, iz) = best_idx;
        end
    end

    error = max(abs(V_grid - V_grid_old), [], 'all'); % compute maximum error between iterations
    iteration = iteration + 1; % increment iteration counter

    if mod(iteration, 50) == 0 
            disp(['Iteration: ', num2str(iteration), ', Error: ', num2str(error)]);
    end
end

time = toc;
disp(['Value function iteration converged in ', num2str(iteration), ' iterations and took ', num2str(time), ' seconds']);

%% Plotting Value and consumption functions
% value function
fig = figure;
theme(fig, "light");
plot(a_grid, V_grid, 'LineWidth', 2);
xlabel('Wealth (a)');
ylabel('Value Function (V)');
title('Value Function Iteration ');
grid on;    
saveas(gcf, 'figs/Value_Function_IFP.png');

% consumption
fig2 = figure;
theme(fig2, "light");                       
c_policy = zeros(Params.n_a, Params.n_z);             % FIX 2: compute consumption policy
for iz = 1:Params.n_z
    z = z_grid(iz);
    for ia = 1:Params.n_a
        a = a_grid(ia);
        a_next = a_grid(policy_grid(ia, iz));         % recover next-period wealth from policy
        c_policy(ia, iz) = (1 + Params.r) * a + z - a_next;
    end
end
plot(a_grid, c_policy, 'LineWidth', 2);               
xlabel('Wealth (a)');                                  
ylabel('Consumption (c)');
title('Consumption Policy Function');
legend(arrayfun(@(i) sprintf('z_%d', i), 1:Params.n_z, 'UniformOutput', false), 'Location', 'best');
grid on;    
saveas(gcf, 'figs/Consumption_Policy_IFP.png');    

% Expected next period wealth against current wealth
fig3 = figure;
theme(fig3, "light");

E_a_next = zeros(Params.n_a, 1);  % initialize expected next period wealth

for ia = 1:Params.n_a
    for iz = 1:Params.n_z
        a_next = a_grid(policy_grid(ia, iz));           % optimal a' for this (a, iz) state
        E_a_next(ia) = E_a_next(ia) + z_grid(iz) * a_next; % weight by stationary prob of z
    end
end

plot(a_grid, E_a_next, 'b-', 'LineWidth', 2); hold on;
plot(a_grid, a_grid, 'k--', 'LineWidth', 1.5);          % 45 degree line
xlabel('Current Wealth (a)');
ylabel('Expected Next Period Wealth E[a'']');
title('Expected Next Period Wealth vs Current Wealth');
legend('E[a''|a]', '45° line', 'Location', 'best');
grid on;
saveas(gcf, 'figs/Expected_Wealth_IFP.png');

% find target wealth for each z state (where a' = a)
fig4 = figure; theme(fig4, "light"); hold on;
for iz = 1:Params.n_z
    a_next_vec = a_grid(policy_grid(:, iz));
    plot(a_grid, a_next_vec - a_grid, 'LineWidth', 2);
end
yline(0, 'k--', 'LineWidth', 1.5);
xlabel('Wealth (a)'); ylabel("a' - a");
title('Target Wealth: Carroll (1997) Buffer Stock');
legend(arrayfun(@(i) sprintf('z_%d', i), 1:Params.n_z, 'UniformOutput', false));
grid on;
saveas(gcf, 'figs/Target_Wealth_IFP.png');

% does wealth grid bind? 
%% Check if Wealth Grid is Binding

% check lower bound -- are any households choosing a' = a_min?
lower_bind = sum(policy_grid(:) == 1);           % count states where a' is at minimum grid point
disp(['Number of (a,z) states at lower bound: ', num2str(lower_bind), ' out of ', num2str(Params.n_a * Params.n_z)]);

% check upper bound -- are any households choosing a' = a_max?
upper_bind = sum(policy_grid(:) == Params.n_a);  % count states where a' is at maximum grid point
disp(['Number of (a,z) states at upper bound: ', num2str(upper_bind), ' out of ', num2str(Params.n_a * Params.n_z)]);

% fraction of states at each bound
disp(['Fraction at lower bound: ', num2str(lower_bind / (Params.n_a * Params.n_z))]);
disp(['Fraction at upper bound: ', num2str(upper_bind / (Params.n_a * Params.n_z))]);


%% Local Functions

% making consumption function 
function c = consumption(a, z, a_grid, Params)
    c_possible= (1 + Params.r) * a + z - a_grid;
    c = max(0, c_possible);
end

% utilty function
function u = utility(c, Params)
    % Vectorized CRRA utility — handles (Na x 1) input
    gamma = Params.gamma;
    u        = -Inf(size(c));       % default: -Inf for c <= 0
    pos      = c > 0;
    u(pos)   = (c(pos).^(1 - gamma)) / (1 - gamma);
end

% cubic spline interpolation for the value function
function V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
    % Replace -Inf with a large negative number for spline stability
    V_grid_clean = V_grid;
    V_grid_clean(~isfinite(V_grid_clean)) = -1e10;
    
    spline_interp = spline(K_grid, V_grid_clean);
    V_interp = ppval(spline_interp, K_query);
end