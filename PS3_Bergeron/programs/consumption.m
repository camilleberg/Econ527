function c = consumption(a, z, a_next, Params)

% c = consumption(a, z, a_next, Params)
%
% Budget constraint: c = (1+r)*a + exp(z) - a'
% z is LOG income so actual income is exp(z)
%
% Inputs:
%   a      - current wealth (scalar or n_a x 1 vector)
%   z      - current log-income (scalar) — will be exponentiated
%   a_next - next-period wealth choice (scalar or n_a x 1 vector)
%   Params - parameter struct

    c = (1 + Params.r) * a + exp(z) - a_next;
end