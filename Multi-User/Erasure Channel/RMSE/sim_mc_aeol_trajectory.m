function [Sim_AEoL_NP, Sim_AEoL_P, Sim_Sys_NP, Sim_Sys_P, Sim_Thr_NP, Sim_Thr_P] = sim_mc_aeol_trajectory(conf)

    num_L = length(conf.lambda_total_vec);
    
    % Matrizes para guardar resultados
    Sim_AEoL_NP = zeros(conf.num_users, num_L);
    Sim_AEoL_P  = zeros(conf.num_users, num_L);
    Sim_Sys_NP  = zeros(1, num_L);
    Sim_Sys_P   = zeros(1, num_L);
    
    % Vetores para Vazão (Throughput)
    Sim_Thr_NP  = zeros(1, num_L);
    Sim_Thr_P   = zeros(1, num_L);

    for j = 1:num_L
        lambda_tot = conf.lambda_total_vec(j);
        lambda_user = lambda_tot / conf.num_users;
        
        acc_np = zeros(conf.num_users, 1);
        acc_p  = zeros(conf.num_users, 1);
        acc_sys_np = 0;
        acc_sys_p  = 0;
        
        acc_count_np = 0;
        acc_count_p  = 0;
        
        for k_mc = 1:conf.N_MC
            % 1. Gerar Chegadas baseadas em TEMPO (não mais N fixo)
            all_arrivals = [];
            
            % Estimatiza segura de quantos pacotes gerar para cobrir Sim_Time
            % Adicionamos margem de segurança de 20% + 50 pacotes
            est_n_pkts = ceil(conf.Sim_Time * lambda_user * 1.2) + 50; 
            
            for u = 1:conf.num_users
                dt = exprnd(1/lambda_user, est_n_pkts, 1); 
                t_arr = cumsum(dt);
                
                % Filtrar apenas o que chegou DENTRO do tempo de simulação
                valid_idx = t_arr <= conf.Sim_Time;
                t_arr = t_arr(valid_idx);
                
                all_arrivals = [all_arrivals; t_arr, ones(length(t_arr), 1)*u];
            end
            [all_arrivals, sort_idx] = sortrows(all_arrivals, 1);
            
            n_pkts = size(all_arrivals, 1);
            
            % Se não houve chegadas (lambda muito baixo ou tempo curto), pular
            if n_pkts == 0
                continue; 
            end
            
            % 2. Gerar Serviços (Erasure Channel)
            services = get_erasure_service(n_pkts, conf.k, conf.delta, conf.Symbol_Duration);
            
            % 3. Gerar Erros de Medição
            meas_errors = normrnd(0, conf.sensor_sigma, n_pkts, 1);
            
            % --- NP Logic ---
            % Retorna [AEoL_User, AEoL_Sys, Num_Processed]
            [aeol_np, sys_aeol_np, count_np] = queue_logic_trajectory(all_arrivals, services, ...
                                          meas_errors, conf.velocities, 'NP', conf.Sim_Time);
            acc_np = acc_np + aeol_np;
            acc_sys_np = acc_sys_np + sys_aeol_np;
            acc_count_np = acc_count_np + count_np;
            
            % --- P Logic ---
            [aeol_p, sys_aeol_p, count_p] = queue_logic_trajectory(all_arrivals, services, ...
                                          meas_errors, conf.velocities, 'P', conf.Sim_Time);
            acc_p = acc_p + aeol_p;
            acc_sys_p = acc_sys_p + sys_aeol_p;
            acc_count_p = acc_count_p + count_p;
        end
        
        % Médias
        Sim_AEoL_NP(:, j) = acc_np / conf.N_MC;
        Sim_AEoL_P(:, j)  = acc_p  / conf.N_MC;
        Sim_Sys_NP(j)     = acc_sys_np / conf.N_MC;
        Sim_Sys_P(j)      = acc_sys_p  / conf.N_MC;
        
        % Throughput = Média de pacotes processados / Tempo de Simulação
        avg_count_np = acc_count_np / conf.N_MC;
        avg_count_p  = acc_count_p  / conf.N_MC;
        
        Sim_Thr_NP(j) = avg_count_np / conf.Sim_Time;
        Sim_Thr_P(j)  = avg_count_p  / conf.Sim_Time;
    end
end

function S = get_erasure_service(n, k, delta, sym_dur)
    if delta >= 1, delta = 0.999; end 
    failures = nbinrnd(k, 1-delta, n, 1);
    S = (k + failures) * sym_dur;
end