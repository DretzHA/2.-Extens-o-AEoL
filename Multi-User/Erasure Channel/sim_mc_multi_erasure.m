function [Sim_NP, Sim_P, Sim_Sys_NP, Sim_Sys_P] = sim_mc_multi_erasure(conf)

    num_L = length(conf.lambda_total_vec);
    
    % Resultados por usuário
    Sim_NP = zeros(conf.num_users, num_L);
    Sim_P  = zeros(conf.num_users, num_L);
    
    % Resultados do Sistema (Novo)
    Sim_Sys_NP = zeros(1, num_L);
    Sim_Sys_P  = zeros(1, num_L);

    for j = 1:num_L
        lambda_tot = conf.lambda_total_vec(j);
        lambda_user = lambda_tot / conf.num_users;
        
        % Acumuladores
        acc_aoi_np = zeros(conf.num_users, 1);
        acc_aoi_p  = zeros(conf.num_users, 1);
        
        acc_sys_np = 0; % Sistema NP
        acc_sys_p  = 0; % Sistema P
        
        for k_mc = 1:conf.N_MC
            % Generate arrivals
            all_arrivals = [];
            for u = 1:conf.num_users
                dt = exprnd(1/lambda_user, conf.N_updates, 1); 
                t_arr = cumsum(dt); 
                all_arrivals = [all_arrivals; t_arr, ones(conf.N_updates, 1)*u];
            end
            
            all_arrivals = sortrows(all_arrivals, 1);
            n_pkts = size(all_arrivals, 1);
            services_needed = get_erasure_service(n_pkts, conf.k, conf.delta);
            services = services_needed*conf.Symbol_Duration;
            
            % --- NP Logic ---
            [aoi_np, sys_np] = queue_logic_multi(all_arrivals, services, conf.num_users, 'NP');
            acc_aoi_np = acc_aoi_np + aoi_np;
            acc_sys_np = acc_sys_np + sys_np;
            
            % --- P Logic ---
            [aoi_p, sys_p] = queue_logic_multi(all_arrivals, services, conf.num_users, 'P');
            acc_aoi_p = acc_aoi_p + aoi_p;
            acc_sys_p = acc_sys_p + sys_p;
        end
        
        % Mean Values (Users)
        mean_aoi_np = acc_aoi_np / conf.N_MC;
        mean_aoi_p  = acc_aoi_p  / conf.N_MC;
        
        % Mean Values (System)
        mean_aoi_Sys_NP = acc_sys_np / conf.N_MC;
        mean_aoi_Sys_P  = acc_sys_p  / conf.N_MC;
        
        % Apply Velocities
        for u = 1:conf.num_users
            v = conf.velocities(u);
            Sim_NP(u, j) = v * mean_aoi_np(u);
            Sim_P(u, j)  = v * mean_aoi_p(u);
        end
        
        Sim_Sys_NP(j) = v * mean_aoi_Sys_NP;
        Sim_Sys_P(j)  = v * mean_aoi_Sys_P ;

        
        % Nota: Para o Sistema, a "velocidade" é ambígua se os usuários tiverem v diferentes.
        % Se todos v=1, ok. Se v diferentes, o "AoI do Sistema" geralmente é puramente temporal (v=1).
        % Se quiser ponderar, precisaria definir uma velocidade média do sistema. 
        % Aqui deixarei puro (tempo).
    end
end

%% FUNCTION ERASURE CHANNEL (Mantida igual)
function S = get_erasure_service(n_packets, k, delta)
    prob_sucess = 1 - delta;
    if prob_sucess >= 1 || delta == 0
        S = ones(n_packets, 1) * k;
    else
        n_failures = nbinrnd(k, prob_sucess, n_packets, 1);
        S = k + n_failures;
    end
end