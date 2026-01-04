function [Sim_NP, Sim_P] = sim_mc_multi(conf)
% SIM_MC_MULTI Gera tráfego de múltiplas fontes e executa a simulação.

    num_L = length(conf.lambda_total_vec);
    Sim_NP = zeros(conf.num_users, num_L);
    Sim_P  = zeros(conf.num_users, num_L);

    for j = 1:num_L
        lambda_tot = conf.lambda_total_vec(j);
        
        % Divide o lambda total pelos utilizadores (Simetria)
        lambda_user = lambda_tot / conf.num_users;
        
        acc_aoi_np = zeros(conf.num_users, 1);
        acc_aoi_p  = zeros(conf.num_users, 1);
        
        for k = 1:conf.N_MC
            % 1. GERAR TRÁFEGO MULTI-FONTE
            all_arrivals = [];
            
            for u = 1:conf.num_users
                % Gera chegadas para o utilizador u
                dt = exprnd(1/lambda_user, conf.N_updates, 1); %tempo entre pacotes
                t_arr = cumsum(dt); %tempo de chegada
                
                % Cria matriz: [Tempo_Chegada, UserID]
                all_arrivals = [all_arrivals; t_arr, ones(conf.N_updates, 1)*u];
            end
            
            % FUNDIR E ORDENAR (Passo Crítico para Multi-User)
            % Ordena todos os pacotes cronologicamente para simular a fila única
            all_arrivals = sortrows(all_arrivals, 1);
            
            % Gera tempos de serviço para o fluxo agregado
            n_pkts = size(all_arrivals, 1);
            services = exprnd(1/conf.mu, n_pkts, 1); %tempo de serviço
            
            % 2. PROCESSAR FILAS (Chama a lógica externa)
            
            % LCFS Non-Preemptive
            aoi_np = queue_logic_multi(all_arrivals, services, conf.num_users, 'NP');
            acc_aoi_np = acc_aoi_np + aoi_np;
            
            % LCFS Preemptive
            aoi_p = queue_logic_multi(all_arrivals, services, conf.num_users, 'P');
            acc_aoi_p = acc_aoi_p + aoi_p;
        end
        
        % Médias e Conversão para AEoL
        mean_aoi_np = acc_aoi_np / conf.N_MC;
        mean_aoi_p  = acc_aoi_p  / conf.N_MC;
        
        for u = 1:conf.num_users
            v = conf.velocities(u);
            Sim_NP(u, j) = v * mean_aoi_np(u);
            Sim_P(u, j)  = v * mean_aoi_p(u);
        end
    end
end