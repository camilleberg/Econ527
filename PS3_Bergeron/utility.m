function u = utility(c, Params)

%  u = utility(c, Params)
%
% CRRA utility function
% c: consumption level (scalar or vector)
% Params: struct of model parameters (e.g.

    gamma = Params.gamma;
    u        = -Inf(size(c));       % default: -Inf for c <= 0
    pos      = c > 0;
    u(pos)   = (c(pos).^(1 - gamma)) / (1 - gamma);
end
