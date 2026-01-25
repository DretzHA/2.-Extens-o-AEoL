function [Sim_AEoL_NP, Sim_AEoL_P, Sim_Sys_NP, Sim_Sys_P, ...
          Sim_Thr_NP, Sim_Thr_P, ...
          Sim_User_Counts_NP, Sim_User_Counts_P] = sim_mc_aeol_trajectory(conf)

    num_L = length(conf.lambda_total_vec); %tamanho do vetor das taxas
    N_u   = conf.num_users; %número de usuários atuais
    
    % Matrizes de Saída do AEoL
    Sim_AEoL_NP = zeros(N_u, num_L); Sim_AEoL_P  = zeros(N_u, num_L); %por ususario
    Sim_Sys_NP  = zeros(1, num_L);   Sim_Sys_P   = zeros(1, num_L); %sistema (junto)
    Sim_Thr_NP  = zeros(1, num_L);   Sim_Thr_P   = zeros(1, num_L); %teórico
    
    % Matrizes de Tempo/Contagem de atualizaçoes
    Sim_Act_Time_NP = zeros(N_u, num_L);    Sim_Act_Time_P  = zeros(N_u, num_L);
    Sim_User_Counts_NP = zeros(N_u, num_L); Sim_User_Counts_P  = zeros(N_u, num_L);

    for j = 1:num_L
        lambda_tot = conf.lambda_total_vec(j); %taxa de atualização atual
        lambda_user = lambda_tot / N_u; %taxa por usuário (dividida igualmente)
        
        % Acumuladores do AEoL
        acc_np = zeros(N_u, 1); acc_p = zeros(N_u, 1); % por user
        acc_sys_np = 0; acc_sys_p = 0; % sistema
        acc_count_total_np = 0; acc_count_total_p = 0; % processados total
        
        acc_act_np = zeros(N_u, 1); acc_act_p = zeros(N_u, 1); % tempo atuação
        acc_usr_cnt_np = zeros(N_u, 1); acc_usr_cnt_p = zeros(N_u, 1); % contagem por user
        
        for k_mc = 1:conf.N_MC
            % Gerar Temps de Chegadas e Serviços
            all_arrivals = [];

            % estima total de pacotes a ser enviado com base no tempo
            % de simulação para alocar tamanho e gerar os arrival
            est_n_pkts = ceil(conf.Sim_Time * lambda_user * 1.2) + 50; 

            for u = 1:N_u
                %tempos de arrival por distribuição exponencial/poisson
                dt = exprnd(1/lambda_user, est_n_pkts, 1); 
                t_arr = cumsum(dt); %vetor de tempo de arrival
                valid_idx = t_arr <= conf.Sim_Time; %verifica se o tempo está abaixo do tempo simulação
                t_arr = t_arr(valid_idx);
                all_arrivals = [all_arrivals; t_arr, ones(length(t_arr), 1)*u]; 
            end

            [all_arrivals, ~] = sortrows(all_arrivals, 1); %sort/organiza pelos tempos
            n_pkts = size(all_arrivals, 1); %verifica quantidade de pacotes
            
            if n_pkts == 0
                continue; 
            end

            %tempo de serviço exponencial para teste de fila M/M/1
            %services = exprnd(1/conf.mu, n_pkts, 1);
            
            %tempo de serviço considerando canal erasure
            services = get_erasure_service(n_pkts, conf.k, conf.delta, conf.Symbol_Duration); 

            %erro RMSE a partir de uma distribuição normal para cada
            %atualização
            meas_errors = normrnd(0, conf.sensor_sigma, n_pkts, 1);
            
            % --- Lógica NP ---
            [aeol_np, sys_aeol_np, count_np, usr_cnt_np] = queue_logic_trajectory(...
                all_arrivals, services, meas_errors, conf.velocities, 'NP', conf.Sim_Time);
            
            acc_np = acc_np + aeol_np;
            acc_sys_np = acc_sys_np + sys_aeol_np;
            acc_count_total_np = acc_count_total_np + count_np;
            acc_usr_cnt_np = acc_usr_cnt_np + usr_cnt_np;
            
            % --- P Logic ---
            [aeol_p, sys_aeol_p, count_p, usr_cnt_p] = queue_logic_trajectory(...
                all_arrivals, services, meas_errors, conf.velocities, 'P', conf.Sim_Time);
            
            acc_p = acc_p + aeol_p;
            acc_sys_p = acc_sys_p + sys_aeol_p;
            acc_count_total_p = acc_count_total_p + count_p;
            acc_usr_cnt_p = acc_usr_cnt_p + usr_cnt_p;
        end
        
        % Médias
        Sim_AEoL_NP(:, j) = acc_np / conf.N_MC;
        Sim_AEoL_P(:, j)  = acc_p  / conf.N_MC;
        Sim_Sys_NP(j)     = acc_sys_np / conf.N_MC;
        Sim_Sys_P(j)      = acc_sys_p  / conf.N_MC;
        
        Sim_Thr_NP(j) = (acc_count_total_np / conf.N_MC) / conf.Sim_Time;
        Sim_Thr_P(j)  = (acc_count_total_p / conf.N_MC)  / conf.Sim_Time;
        
        % Métricas por Usuário
        Sim_Act_Time_NP(:, j) = acc_act_np / conf.N_MC;
        Sim_Act_Time_P(:, j)  = acc_act_p  / conf.N_MC;
        
        Sim_User_Counts_NP(:, j) = acc_usr_cnt_np / conf.N_MC;
        Sim_User_Counts_P(:, j)  = acc_usr_cnt_p  / conf.N_MC;
    end
end

%Função para obter tempo de serviço de canal erasure
function S = get_erasure_service(n, k, delta, sym_dur)

    if delta >= 1
       delta = 0.999; 
    end 

    failures = nbinrnd(k, 1-delta, n, 1); %quantidade de retransmissões são um binomial negativo - Yates 2017
    S = (k + failures) * sym_dur; %tempo de serviço total
end