function c = consumption(a, z, a_grid, Params)

% c = consumption(a, z, a_grid, Params)
% calculates consumption for given current wealth (a), income (z), and next period wealth choices (a
%
% a: current wealth
% z: current income
% a_grid: grid of next period wealth choices (a')
% Params: struct of model parameters (e.g. beta, signa, etc.)

    c_possible = (1 + Params.r) * a + z - a_grid;
    c = max(0, c_possible);
end
