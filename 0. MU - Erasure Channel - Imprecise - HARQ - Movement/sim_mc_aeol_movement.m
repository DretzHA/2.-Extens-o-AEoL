function [Sim_AEoL_NP, Sim_AEoL_P, Sim_Sys_NP, Sim_Sys_P, ...
          Sim_User_Counts_NP, Sim_User_Counts_P] = sim_mc_aeol_movement(conf)

    num_L = length(conf.lambda_total_vec);
    N_u   = conf.num_users;
    
    % Matrizes de Saída
    Sim_AEoL_NP = zeros(N_u, num_L); Sim_AEoL_P  = zeros(N_u, num_L);
    Sim_Sys_NP  = zeros(1, num_L);   Sim_Sys_P   = zeros(1, num_L);
    Sim_User_Counts_NP = zeros(N_u, num_L); Sim_User_Counts_P  = zeros(N_u, num_L);

    for j = 1:num_L
        lambda_tot = conf.lambda_total_vec(j);
        lambda_user = lambda_tot / N_u;
        
        acc_sys_np = 0; acc_sys_p = 0;
        acc_usr_np = zeros(N_u, 1); acc_usr_p = zeros(N_u, 1);
        acc_cnt_np = zeros(N_u, 1); acc_cnt_p = zeros(N_u, 1);

        parfor k = 1:conf.N_MC % 
            
            % 1. Gerar Trajetória RWP 2D (Comum para P e NP nesta iteração)
            % Retorna [t_start, t_end, start_x, start_y, vel_x, vel_y]
            path_segments = generate_rwp_path(conf.Sim_Time, conf.area_side, conf.velocity);
            
            % 2. Gerar Chegadas e Erros de Medição
            all_arrivals = [];
            all_meas_errors = [];
            
            for u = 1:N_u
                % Chegadas Poisson
                dt = exprnd(1/lambda_user, ceil(conf.Sim_Time * lambda_user * 2) + 50, 1);
                t_arr = cumsum(dt);
                t_arr = t_arr(t_arr <= conf.Sim_Time);
                
                % Erro de Medição (x, y) ~ N(0, sigma)
                err_u = conf.sensor_sigma * randn(length(t_arr), 2);
                
                all_arrivals = [all_arrivals; t_arr, ones(length(t_arr), 1)*u];
                all_meas_errors = [all_meas_errors; err_u];
            end
            
            % Ordenar chegadas
            [~, idx_sort] = sort(all_arrivals(:, 1));
            all_arrivals = all_arrivals(idx_sort, :);
            all_meas_errors = all_meas_errors(idx_sort, :);
            
            n_pkts = size(all_arrivals, 1);
            
            % 3. Gerar Tempos de Serviço (Canal Erasure / HARQ)
            transmissions = geornd(1 - conf.delta, n_pkts, 1) + 1; %tentativas até conseguir (pacote inteiro)
            services = transmissions * conf.k * conf.Symbol_Duration;
            
            % --- Simulação LCFS NON-PREEMPTIVE ---
            [aeol_np, sys_aeol_np, ~, usr_cnt_np] = queue_logic_movement(...
                all_arrivals, services, all_meas_errors, path_segments, 'NP', conf.Sim_Time);
            
            acc_usr_np = acc_usr_np + aeol_np;
            acc_sys_np = acc_sys_np + sys_aeol_np;
            acc_cnt_np = acc_cnt_np + usr_cnt_np;
            
            % --- Simulação LCFS PREEMPTIVE ---
            [aeol_p, sys_aeol_p, ~, usr_cnt_p] = queue_logic_movement( ...
                all_arrivals, services, all_meas_errors, path_segments, 'P', conf.Sim_Time);
            
            acc_usr_p = acc_usr_p + aeol_p;
            acc_sys_p = acc_sys_p + sys_aeol_p;
            acc_cnt_p = acc_cnt_p + usr_cnt_p;
        end
        
        % Médias Monte Carlo
        Sim_AEoL_NP(:, j) = acc_usr_np / conf.N_MC;
        Sim_AEoL_P(:, j)  = acc_usr_p  / conf.N_MC;
        Sim_Sys_NP(j)     = acc_sys_np / conf.N_MC;
        Sim_Sys_P(j)      = acc_sys_p  / conf.N_MC;
        
        Sim_User_Counts_NP(:, j) = acc_cnt_np / conf.N_MC;
        Sim_User_Counts_P(:, j)  = acc_cnt_p  / conf.N_MC;
    end
end

function segments = generate_rwp_path(T_max, L, v)
    % Gera segmentos de movimento retilíneo (RWP)
    % segments: [t_start, t_end, x0, y0, vx, vy]
    
    t_curr = 0;
    curr_pos = rand(1, 2) * L; % Posição inicial aleatória
    segments = [];
    
    while t_curr < T_max
        target_pos = rand(1, 2) * L;
        dist = norm(target_pos - curr_pos);
        
        if dist == 0
            continue;
        end
        
        travel_time = dist / v;
        velocity_vec = (target_pos - curr_pos) / travel_time;
        
        t_next = t_curr + travel_time;
        
        % [t_start, t_end, x0, y0, vx, vy]
        segments = [segments; t_curr, t_next, curr_pos(1), curr_pos(2), velocity_vec(1), velocity_vec(2)];
        
        t_curr = t_next;
        curr_pos = target_pos;
    end
end
