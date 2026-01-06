function [Sim_AEoL_NP, Sim_AEoL_P, Sim_Sys_NP, Sim_Sys_P] = sim_mc_aeol_trajectory(conf)

    num_L = length(conf.lambda_total_vec);
    
    Sim_AEoL_NP = zeros(conf.num_users, num_L);
    Sim_AEoL_P  = zeros(conf.num_users, num_L);
    Sim_Sys_NP  = zeros(1, num_L);
    Sim_Sys_P   = zeros(1, num_L);

    for j = 1:num_L
        lambda_tot = conf.lambda_total_vec(j);
        lambda_user = lambda_tot / conf.num_users;
        
        acc_np = zeros(conf.num_users, 1);
        acc_p  = zeros(conf.num_users, 1);
        acc_sys_np = 0;
        acc_sys_p  = 0;
        
        for k_mc = 1:conf.N_MC
            % 1. Gerar Chegadas
            all_arrivals = [];
            for u = 1:conf.num_users
                dt = exprnd(1/lambda_user, conf.N_updates, 1); 
                t_arr = cumsum(dt); 
                all_arrivals = [all_arrivals; t_arr, ones(conf.N_updates, 1)*u];
            end
            [all_arrivals, sort_idx] = sortrows(all_arrivals, 1);
            
            % 2. Gerar Serviços (Erasure Channel)
            n_pkts = size(all_arrivals, 1);
            services = get_erasure_service(n_pkts, conf.k, conf.delta, conf.Symbol_Duration);
            
            % 3. Gerar Erros de Medição (Noise Samples)
            % Cada pacote carrega uma medição que é: Pos_Real + Erro
            % O erro vem de uma Normal(0, sigma). Pode ser positivo ou negativo.
            % Se quiser simular RMSE variável por pacote, altere o sigma aqui dentro.
            meas_errors = normrnd(0, conf.sensor_sigma, n_pkts, 1);
            
            % --- NP Logic ---
            [aeol_np, sys_aeol_np] = queue_logic_trajectory(all_arrivals, services, ...
                                          meas_errors, conf.velocities, 'NP');
            acc_np = acc_np + aeol_np;
            acc_sys_np = acc_sys_np + sys_aeol_np;
            
            % --- P Logic ---
            [aeol_p, sys_aeol_p] = queue_logic_trajectory(all_arrivals, services, ...
                                          meas_errors, conf.velocities, 'P');
            acc_p = acc_p + aeol_p;
            acc_sys_p = acc_sys_p + sys_aeol_p;
        end
        
        % Médias
        Sim_AEoL_NP(:, j) = acc_np / conf.N_MC;
        Sim_AEoL_P(:, j)  = acc_p  / conf.N_MC;
        Sim_Sys_NP(j)     = acc_sys_np / conf.N_MC;
        Sim_Sys_P(j)      = acc_sys_p  / conf.N_MC;
    end
end

function S = get_erasure_service(n, k, delta, sym_dur)
    % Gera tempo de serviço considerando retransmissões geométricas (Erasure)
    if delta >= 1, delta = 0.999; end 
    failures = nbinrnd(k, 1-delta, n, 1);
    S = (k + failures) * sym_dur;
end