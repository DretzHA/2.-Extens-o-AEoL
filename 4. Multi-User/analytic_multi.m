function [Res_NP, Res_P] = analytic_multi(conf)
% ANALYTIC_MULTI Calcula o AEoL teórico baseado em Yates & Kaul (2019).
% Saída: Matrizes [Num_Users x Num_Lambdas]

    num_L = length(conf.lambda_total_vec);
    Res_NP = zeros(conf.num_users, num_L);
    Res_P  = zeros(conf.num_users, num_L);

    for j = 1:num_L
        lambda_tot = conf.lambda_total_vec(j);
        rho = lambda_tot / conf.mu;
        
        % Assumindo divisão igual da taxa de chegada (Simétrico)
        lambda_i = lambda_tot / conf.num_users;
        
        % Taxa de utilização individual (necessária para o Teorema 5)
        rho_i = lambda_i / conf.mu;
        
        for u = 1:conf.num_users
            v = conf.velocities(u);
            
            % --- 1. LCFS-W ---
            
            alfa_w_nom = (1+rho+rho^2)^2+2*rho^3;
            alfa_w_denom = (1+rho+rho^2)*(1+rho)^2;
            alfa_w = alfa_w_nom/alfa_w_denom;

            term1 = alfa_w + (1+(rho^2/(1+rho)))/rho_i;

            AoI_NP = term1/conf.mu;
            
            Res_NP(u, j) = v * AoI_NP;
            
            % --- 2. LCFS-S ---

            AoI_P = (1/conf.mu)*(1+rho)*(1/rho_i);
            
            Res_P(u, j) = v * AoI_P;
        end
    end
end