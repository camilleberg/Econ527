function c = u_prime_inv(mu, Params)

% c = u_prime_inv(mu, Params)
%
% Returns the inverse of the marginal utility function u'(c) = c^{-gamma}.
%
% Since  u'(c) = c^{-gamma}  =>  c = mu^{-1/gamma}
%
% Inputs:
%   mu     - marginal utility value(s) (scalar or array), must be > 0
%   Params - (optional) struct with field .gamma
%
% Output:
%   c      - consumption level(s) satisfying u'(c) = mu
    c = mu .^ (-1 / Params.gamma);
end