function u_prime_inv = utility_prime_inv(u, Params) % (Na x 1)  

% u_prime_inv = utility_primeinv(u, Params)
%
% calcualtes inverse of CRRA utility
% 
% u: utility level
% Params: struc of model

    u_prime_inv = u.^(-1/Params.gamma); % inverse of CRRA utility
end
