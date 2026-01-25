function [AEoL_woP, AEoL_withP] = compute_analytical_results(velocities, mu, rho_vec, lambda_vec)

    n_v = length(velocities);
    n_l = length(lambda_vec);
    AEoL_woP = zeros(n_v, n_l);
    AEoL_withP = zeros(n_v, n_l);

    for j = 1:n_l
        lambda = lambda_vec(j);
        cur_rho = rho_vec(j);
        
        for i = 1:n_v
            v = velocities(i);
            
            % --- 1. M/M/1 LCFS Sem Preempção ---
            num = cur_rho^3*(1+cur_rho)^2 + cur_rho^4 + (1+cur_rho)^3;
            den = cur_rho^4*(1+cur_rho) + cur_rho*(1+cur_rho)^3;
            termo_woP = 1 + num/den;
            AEoL_woP(i, j) = (v / mu) * termo_woP;
            
            % --- 2. M/M/1 LCFS Com Preempção ---
            AEoL_withP(i, j) = v * (1/lambda + 1/mu);
        end
    end
end