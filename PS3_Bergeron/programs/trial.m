%% VFI with Cubic Interpolation: Income Fluctuation Problem
clear; clc;

% 1. Parameters
beta  = 0.95;       % Discount factor
gamma = 2.0;        % Risk aversion (CRRA)
R     = 1.02;       % Gross interest rate
y_states = [0.8, 1.2]; % Income states (Low, High)
Pi = [0.9 0.1; 0.1 0.9]; % Transition matrix

% 2. Asset Grid
a_min = 0; a_max = 20; n_a = 50;
a_grid = linspace(a_min, a_max, n_a)';

% 3. Initialize Value Function
V = zeros(n_a, length(y_states));
V_new = V;
policy_a = zeros(n_a, length(y_states));
tol = 1e-6; max_iter = 500; diff = 1; iter = 0;

% 4. VFI Loop
while diff > tol && iter < max_iter
    % Pre-calculate Expected Value for each asset choice
    % This speeds up the inner loop significantly
    EV = V * Pi'; 
    
    for j = 1:length(y_states) % Loop over income
        for i = 1:n_a          % Loop over current assets
            a_curr = a_grid(i);
            y_curr = y_states(j);
            
            % RHS of Bellman: -[u(c) + beta * interp(V)]
            % We minimize the negative to use fminbnd
            obj = @(a_next) -( ((R*a_curr + y_curr - a_next)^(1-gamma)-1)/(1-gamma) ...
                + beta * spline(a_grid, EV(:,j), a_next) );
            
            % Constraints: Cannot borrow past a_min, cannot consume more than wealth
            lb = a_min;
            ub = min(R*a_curr + y_curr - 0.01, a_max);
            
            % Continuous optimization for a_prime
            [a_prime, fval] = fminbnd(obj, lb, ub);
            
            V_new(i,j) = -fval;
            policy_a(i,j) = a_prime;
        end
    end
    
    diff = max(abs(V_new(:) - V(:)));
    V = V_new;
    iter = iter + 1;
    fprintf('Iteration %d, Diff: %e\n', iter, diff);
end

% 5. Plot Results
figure;
plot(a_grid, policy_a);
title('Optimal Asset Policy Function (Cubic Interpolation)');
legend('Low Income', 'High Income');
xlabel('Current Assets'); ylabel('Next Period Assets');
