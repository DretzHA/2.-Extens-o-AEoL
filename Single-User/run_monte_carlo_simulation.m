function [AEoL_woP, AEoL_withP] = run_monte_carlo_simulation(velocities, mu, lambda_vec, N_MC, N_updates)
% RUN_MONTE_CARLO_SIMULATION Executa o loop principal de simulação numérica.
% Gera chegadas e serviços aleatórios e chama a função de fila.

    n_v = length(velocities);
    n_l = length(lambda_vec);
    
    AEoL_woP = zeros(n_v, n_l);
    AEoL_withP = zeros(n_v, n_l);

    % Loop pelas Taxas de Chegada
    for j = 1:n_l
        lambda = lambda_vec(j);
        
        avg_AoI_woP_acc = 0;
        avg_AoI_withP_acc = 0;
        
        % Loop Monte Carlo (Realizações)
        for k = 1:N_MC
            % Gera tempos aleatórios (M/M/1)
            inter_T = exprnd(1/lambda, N_updates, 1);
            Serv_T  = exprnd(1/mu, N_updates, 1);
            
            % Chama a função da fila (arquivo separado)
            
            % 1. Sem Preempção
            avg_AoI_woP_acc = avg_AoI_woP_acc + ...
                calculate_aoi_lcfs(inter_T, Serv_T, 'NonPreemptive');
                
            % 2. Com Preempção
            avg_AoI_withP_acc = avg_AoI_withP_acc + ...
                calculate_aoi_lcfs(inter_T, Serv_T, 'Preemptive');
        end
        
        % Médias de AoI
        mean_AoI_woP = avg_AoI_woP_acc / N_MC;
        mean_AoI_withP = avg_AoI_withP_acc / N_MC;
        
        % Cálculo do AEoL Simulado (AEoL = v * AoI)
        for i = 1:n_v
            v = velocities(i);
            AEoL_woP(i, j) = v * mean_AoI_woP;
            AEoL_withP(i, j) = v * mean_AoI_withP;
        end
    end
end