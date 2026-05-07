function c = utility_prime_inv(u, Params) % (Na x 1)  

% u_prime_inv = utility_primeinv(u, Params)
%
% calcualtes inverse of CRRA utility
% 
% u: utility level
% Params: struc of model
% since u' = c^(-gamma), u'^-1 = c = u^(-1/gamma)

    c= u.^(-1/Params.gamma); % inverse of CRRA utility
end
