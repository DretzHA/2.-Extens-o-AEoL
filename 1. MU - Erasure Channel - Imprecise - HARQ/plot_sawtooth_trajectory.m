function plot_sawtooth_trajectory(conf, n_users_to_plot)
    % Gera uma realização única e plota o AGE OF INFORMATION (AoI).
    % Compara:
    % 1. LCFS Non-Preemptive (NP)
    % 2. LCFS Preemptive (P)
    % Considera canais com Erasure (HARQ/Retransmissão via tempo de serviço estendido).
    
    fprintf('Gerando gráfico Sawtooth (AoI - NP vs P)...\n');
    
    % --- 1. Gerar dados para UMA realização ---
    idx_lambda = 7; 
    lambda_tot = conf.lambda_total_vec(idx_lambda); 
    lambda_user = lambda_tot / conf.num_users;
    
    all_arrivals = [];
    % Gera pacotes suficientes para cobrir o tempo de simulação
    est_n_pkts = ceil(conf.Sim_Time * lambda_user * 2) + 50; 
    
    for u = 1:conf.num_users
        dt = exprnd(1/lambda_user, est_n_pkts, 1);
        t_arr = cumsum(dt);
        valid = t_arr <= conf.Sim_Time;
        t_arr = t_arr(valid);
        all_arrivals = [all_arrivals; t_arr, ones(length(t_arr), 1)*u];
    end
    [all_arrivals, ~] = sortrows(all_arrivals, 1);
    
    n_pkts = size(all_arrivals, 1);
    if n_pkts == 0, warning('Nenhum pacote gerado.'); return; end
    
    % Serviços (com Erasure/HARQ) e Erros (não usados no AoI, mas mantidos na estrutura)
    services = get_erasure_service(n_pkts, conf.k, conf.delta, conf.Symbol_Duration);
    meas_errors = normrnd(0, conf.sensor_sigma, n_pkts, 1);
    
    % Dados Completos [Chegada, Usuario, Servico, ErroMedicao]
    data_full = [all_arrivals, services, meas_errors];
    
    % --- 2. Processar Filas ---
    completed_np = run_lcfs_np_plot(data_full, conf.Sim_Time);
    completed_p  = run_lcfs_p_plot(data_full, conf.Sim_Time);
    
    if n_users_to_plot > conf.num_users, n_users_to_plot = conf.num_users; end
    
    % =====================================================================
    % PLOTAGEM
    % =====================================================================
    figure('Name', 'AoI Sawtooth: NP vs P', 'Color', 'w', 'Position', [100, 100, 1200, 800]);
    
    % Limites de Zoom para visualização (ex: primeiros 50s ou tudo)
    time_limit = min(50, conf.Sim_Time);
    
    % --- PLOT 1: LCFS - Não Preemptivo (NP) ---
    subplot(2, 1, 1); hold on; grid on; box on;
    title('Age of Information (AoI) - LCFS Non-Preemptive', 'Interpreter', 'latex', 'FontSize', 14);
    
    [t_vec_np, aoi_np] = generate_aoi_curve(completed_np, conf.Sim_Time);
    
    % Plot da curva de AoI
    plot(t_vec_np, aoi_np, '-', 'Color', 'b', 'LineWidth', 1.2, 'DisplayName', 'AoI System (NP)');
    
    % Marcadores de entrega (Updates)
    if ~isempty(completed_np)
        % No eixo X é o tempo de entrega (FinishTime), no Y é o AoI naquele instante (Finish - Gen)
        finish_times = completed_np(:,1);
        gen_times = completed_np(:,2);
        aoi_peaks = finish_times - gen_times; 
        plot(finish_times, aoi_peaks, 'v', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', ...
             'MarkerSize', 5, 'DisplayName', 'Updates');
    end
    
    ylabel('Age of Information (s)', 'Interpreter', 'latex');
    xlim([0, time_limit]);
    ylim([0, max(aoi_np)*1.05]);
    legend('Location', 'best');

    % --- PLOT 2: LCFS - Preemptivo (P) ---
    subplot(2, 1, 2); hold on; grid on; box on;
    title('Age of Information (AoI) - LCFS Preemptive', 'Interpreter', 'latex', 'FontSize', 14);
    
    [t_vec_p, aoi_p] = generate_aoi_curve(completed_p, conf.Sim_Time);
    
    % Plot da curva de AoI
    plot(t_vec_p, aoi_p, '-', 'Color', 'r', 'LineWidth', 1.2, 'DisplayName', 'AoI System (P)');
    
    % Marcadores de entrega
    if ~isempty(completed_p)
        finish_times = completed_p(:,1);
        gen_times = completed_p(:,2);
        aoi_peaks = finish_times - gen_times; 
        plot(finish_times, aoi_peaks, 'v', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', ...
             'MarkerSize', 5, 'DisplayName', 'Updates');
    end
    
    xlabel('Time (s)', 'Interpreter', 'latex');
    ylabel('Age of Information (s)', 'Interpreter', 'latex');
    xlim([0, time_limit]);
    % Ajusta escala Y baseada no NP para facilitar comparação visual, ou usa auto
    % ylim([0, max(aoi_np)*1.05]); 
    legend('Location', 'best');
    
end

% =========================================================================
%  HELPERS
% =========================================================================

function [t_out, aoi_out] = generate_aoi_curve(completed_data, T_max)
    % Gera vetores de tempo e AoI
    % AoI(t) = t - t_generation_of_last_received_packet
    
    % Se nada foi entregue, assume geração t=0
    if isempty(completed_data)
        t_out = [0, T_max];
        aoi_out = [0, T_max];
        return;
    end

    data = sortrows(completed_data, 1); % Ordenar por tempo de entrega (Finish Time)
    
    t_out = [];
    aoi_out = [];
    
    current_t = 0;
    last_gen = 0; % Assumimos que no t=0 temos informação fresca (AoI=0) ou gen=0
    
    dt_res = 0.05; % Resolução fina
    
    finish_times = data(:, 1);
    gen_times    = data(:, 2);
    
    for k = 1:length(finish_times)
        next_t = finish_times(k);
        if next_t > T_max, next_t = T_max; end
        
        if next_t > current_t
            % Segmento entre a atualização anterior e a atual
            t_seg = current_t : dt_res : next_t;
            if t_seg(end) ~= next_t, t_seg = [t_seg, next_t]; end
            
            % Cálculo do AoI: Crescimento Linear com inclinação 1
            seg_aoi = t_seg - last_gen;
            
            t_out = [t_out, t_seg];
            aoi_out = [aoi_out, seg_aoi];
        end
        
        % Update State
        current_t = next_t;
        last_gen = gen_times(k); % O sistema agora conhece o pacote gerado em gen_times(k)
        
        if current_t >= T_max, break; end
    end
    
    % Preencher o restante do tempo até T_max
    if current_t < T_max
        t_seg = current_t : dt_res : T_max;
        seg_aoi = t_seg - last_gen;
        t_out = [t_out, t_seg];
        aoi_out = [aoi_out, seg_aoi];
    end
end

function S = get_erasure_service(n, k, delta, sym_dur)
    % Gera tempo de serviço baseado em Binomial Negativa (Erasure Channel)
    if delta >= 1, delta = 0.999; end 
    failures = nbinrnd(k, 1-delta, n, 1);
    S = (k + failures) * sym_dur;
end

function completed = run_lcfs_np_plot(data, T_max)
    % data: [Arrival, User, Service, Error]
    N = size(data, 1);
    completed = [];
    served_mask = false(N, 1);
    time_now = 0;
    
    while true
        if time_now >= T_max, break; end
        % Encontra pacotes que já chegaram e não foram servidos
        queue = find(data(:,1) <= time_now & ~served_mask);
        
        if isempty(queue)
            % Ocioso: avança para a próxima chegada
            upcoming = find(data(:,1) > time_now & ~served_mask, 1);
            if isempty(upcoming), break; end
            time_now = data(upcoming, 1); continue;
        end
        
        % LCFS: Pega o último que chegou
        target = queue(end);
        
        % Drop logic (NP): Todos os outros na fila são descartados se for buffer=1
        % Se buffer infinito mas política LCFS pura, eles ficam lá, mas 
        % geralmente em AoI assume-se buffer de tamanho 1 ou substituição.
        % Aqui assumimos fila com "discard old" implícito ao pegar o queue(end)
        % e marcar os outros como servidos/dropados?
        % No código original: "if length(queue) > 1, served_mask(...) = true;"
        % Isso implica DROP dos pacotes velhos (Gerenciamento de fila LCFS com descarte).
        if length(queue) > 1
            served_mask(queue(1:end-1)) = true; 
        end
        
        fin_t = time_now + data(target, 3);
        
        if fin_t <= T_max
            % [Finish, Gen, User, Error]
            completed = [completed; fin_t, data(target, 1), data(target, 2), data(target, 4)];
            time_now = fin_t;
        else
            time_now = fin_t; 
        end
        served_mask(target) = true;
    end
end

function completed = run_lcfs_p_plot(data, T_max)
    % Simulação LCFS Preemptiva
    % Se chegar um pacote novo enquanto serve, o atual é descartado.
    
    N = size(data, 1);
    completed = [];
    
    % Iterar sobre todos os pacotes
    for k = 1:N-1
        arr_curr = data(k, 1);
        arr_next = data(k+1, 1);
        serv_curr = data(k, 3);
        
        finish_time = arr_curr + serv_curr;
        
        % Condição de Sobrevivência na Preempção:
        % O pacote deve terminar ANTES que o próximo chegue.
        if finish_time <= arr_next && finish_time <= T_max
             % [Finish, Gen, User, Error]
             completed = [completed; finish_time, arr_curr, data(k, 2), data(k, 4)];
        else
            % Foi preemptado pelo pacote k+1 (Drop)
        end
    end
    
    % Tratar o último pacote (ninguém para preemptar ele)
    last_idx = N;
    finish_time = data(last_idx, 1) + data(last_idx, 3);
    if finish_time <= T_max
        completed = [completed; finish_time, data(last_idx, 1), data(last_idx, 2), data(last_idx, 4)];
    end
end