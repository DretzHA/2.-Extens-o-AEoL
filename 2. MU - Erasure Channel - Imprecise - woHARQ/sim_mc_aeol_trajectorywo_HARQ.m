function [Sim_AEoL_NP, Sim_AEoL_P, Sim_Sys_NP, Sim_Sys_P, ...
          Sim_Thr_NP, Sim_Thr_P, ...
          Sim_Act_Time_NP, Sim_Act_Time_P, ...
          Sim_User_Counts_NP, Sim_User_Counts_P] = sim_mc_aeol_trajectorywo_HARQ(conf)

    num_L = length(conf.lambda_total_vec);
    N_u   = conf.num_users;
    
    % Matrizes de Saída
    Sim_AEoL_NP = zeros(N_u, num_L); Sim_AEoL_P  = zeros(N_u, num_L);
    Sim_Sys_NP  = zeros(1, num_L);   Sim_Sys_P   = zeros(1, num_L);
    Sim_Thr_NP  = zeros(1, num_L);   Sim_Thr_P   = zeros(1, num_L);
    
    Sim_Act_Time_NP = zeros(N_u, num_L);
    Sim_Act_Time_P  = zeros(N_u, num_L);
    Sim_User_Counts_NP = zeros(N_u, num_L);
    Sim_User_Counts_P  = zeros(N_u, num_L);

    for j = 1:num_L
        lambda_tot = conf.lambda_total_vec(j);
        lambda_user = lambda_tot / N_u;
        
        acc_np = zeros(N_u, 1); acc_p = zeros(N_u, 1);
        acc_sys_np = 0;         acc_sys_p = 0;
        
        acc_count_total_np = 0; acc_count_total_p = 0;
        acc_act_np = zeros(N_u, 1); acc_act_p = zeros(N_u, 1);
        acc_usr_cnt_np = zeros(N_u, 1); acc_usr_cnt_p = zeros(N_u, 1);

        for k_mc = 1:conf.N_MC
            % --- 1. Gerar Chegadas (Poisson) ---
            all_arrivals = [];
            % Estimar numero de pacotes para alocar memória
            est_n_pkts = ceil(conf.Sim_Time * lambda_user * 1.5) + 50; 
            
            for u = 1:N_u
                dt = exprnd(1/lambda_user, est_n_pkts, 1);
                t_arr = cumsum(dt);
                valid = t_arr <= conf.Sim_Time;
                t_arr = t_arr(valid);
                % Col 1: Tempo, Col 2: UserID
                all_arrivals = [all_arrivals; t_arr, ones(length(t_arr), 1)*u];
            end
            % Ordenar chegadas globalmente
            [all_arrivals, ~] = sortrows(all_arrivals, 1);
            n_pkts = size(all_arrivals, 1);
            
            if n_pkts == 0, continue; end
            
            % --- 2. Gerar Serviço e Sucesso (SEM HARQ) ---
            % Serviço Fixo: K * Duração
            fixed_service = conf.k * conf.Symbol_Duration;
            services = ones(n_pkts, 1) * fixed_service;
            
            % Canal Erasure
            pkt_success_prob = (1 - conf.delta);       
            is_successful = rand(n_pkts, 1) < pkt_success_prob;
            
            % --- 3. Erro de Medição ---
            meas_errors = normrnd(0, conf.sensor_sigma, n_pkts, 1);
            
            % --- 4. Executar Lógica de Fila ---
            
            % NP Logic
            [aeol_np, sys_aeol_np, count_np, act_np, usr_cnt_np] = queue_logic_woHARQ(...
                all_arrivals, services, meas_errors, conf.velocities, 'NP', conf.Sim_Time, is_successful);
            
            acc_np = acc_np + aeol_np;
            acc_sys_np = acc_sys_np + sys_aeol_np;
            acc_count_total_np = acc_count_total_np + count_np;
            acc_act_np = acc_act_np + act_np;
            acc_usr_cnt_np = acc_usr_cnt_np + usr_cnt_np;
            
            % P Logic
            [aeol_p, sys_aeol_p, count_p, act_p, usr_cnt_p] = queue_logic_woHARQ(...
                all_arrivals, services, meas_errors, conf.velocities, 'P', conf.Sim_Time, is_successful);
            
            acc_p = acc_p + aeol_p;
            acc_sys_p = acc_sys_p + sys_aeol_p;
            acc_count_total_p = acc_count_total_p + count_p;
            acc_act_p = acc_act_p + act_p;
            acc_usr_cnt_p = acc_usr_cnt_p + usr_cnt_p;
        end
        
        % Médias Monte Carlo
        Sim_AEoL_NP(:, j) = acc_np / conf.N_MC;
        Sim_AEoL_P(:, j)  = acc_p  / conf.N_MC;
        Sim_Sys_NP(j)     = acc_sys_np / conf.N_MC;
        Sim_Sys_P(j)      = acc_sys_p  / conf.N_MC;
        
        Sim_Thr_NP(j) = (acc_count_total_np / conf.N_MC) / conf.Sim_Time;
        Sim_Thr_P(j)  = (acc_count_total_p / conf.N_MC) / conf.Sim_Time;
        
        Sim_Act_Time_NP(:, j) = acc_act_np / conf.N_MC;
        Sim_Act_Time_P(:, j)  = acc_act_p / conf.N_MC;
        
        Sim_User_Counts_NP(:, j) = acc_usr_cnt_np / conf.N_MC;
        Sim_User_Counts_P(:, j)  = acc_usr_cnt_p / conf.N_MC;
        
    end
end
