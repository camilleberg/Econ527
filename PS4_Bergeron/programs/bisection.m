function r = bisection(r_low, r_high)
% r = bisection(r_low, r_high)
%
% finds equlbrium interest rate

% verifies that the functions are positive 
nad_high = f_nad(z_grid, z_prob, @u_fxn, @u_prime_inv, Params, a_next_grid, a_fine_grid);


    while (r_high - r_low) > 1e-4
        r_mid = (r_low + r_high) / 2; 
        Params.r = r_mid; 
        tic;
        net_assets = f_nad(z_grid, z_prob, @u_fxn, @u_prime_inv, Params, a_next_grid, a_fine_grid);
        time_nad = toc; 
        fprintf('Net assets are %.4f at real interest rate %.4f\n and took %.4f to run \n', net_assets, Params.r, time_nad);
        
        if net_assets > 0
            r_low = r_mid;    % excess saving → r too low, raise the floor
        else
            r_high = r_mid;   % excess borrowing → r too high, lower the ceiling
        end
        r = (r_low + r_high)/2;
    fprintf('\nEquilibrium interest rate: r* = %.6f\n', (r_low + r_high)/2);
    end
end