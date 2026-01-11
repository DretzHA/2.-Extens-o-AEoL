function plot_sawtooth_trajectory(conf, n_users_to_plot)
    % Gera uma realização única e plota:
    % 1. Curva Real (AEoL): Inclui erro de medição (sensor_sigma).
    % 2. Curva Ideal (AoI Scaled): Supõe erro de medição = 0.
    
    fprintf('Gerando gráfico Sawtooth (Real vs Ideal)...\n');
    
    % --- 1. Gerar dados para UMA realização ---
    lambda_tot = conf.lambda_total_vec(end); 
    lambda_user = lambda_tot / conf.num_users;
    
    all_arrivals = [];
    est_n_pkts = ceil(conf.Sim_Time * lambda_user * 1.5) + 50; 
    
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
    
    % Serviços e Erros
    services = get_erasure_service(n_pkts, conf.k, conf.delta, conf.Symbol_Duration);
    meas_errors = normrnd(0, conf.sensor_sigma, n_pkts, 1);
    velocities = conf.velocities;
    
    % Dados Completos
    data_full = [all_arrivals, services, meas_errors];
    
    % Rodar Lógica (Ex: NP)
    completed = run_lcfs_np_plot(data_full, conf.Sim_Time);
    
    if n_users_to_plot > conf.num_users, n_users_to_plot = conf.num_users; end
    
    % =====================================================================
    % PLOTAGEM
    % =====================================================================
    figure('Name', 'Trajectory Error: Real vs Ideal', 'Color', 'w', 'Position', [100, 100, 1000, 700]);
    
    % --- PLOT 1: Visão do Sistema (Real vs Ideal) ---
    subplot(2, 1, 1); hold on; grid on; box on;
    
    % Obter as curvas
    v_mean = mean(velocities);
    [t_vec, err_real, err_ideal] = generate_error_curve(completed, v_mean, conf.Sim_Time);
    
    % 1. Plot Ideal (Fundo, tracejado verde)
    plot(t_vec, err_ideal, '--', 'Color', [0 0.6 0], 'LineWidth', 1.5, ...
         'DisplayName', 'Ideal (No Meas. Error)');
    
    % 2. Plot Real (Frente, sólido preto/azul)
    % Preenchimento opcional para destacar a diferença
    % fill([t_vec, fliplr(t_vec)], [err_real, fliplr(err_ideal)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    plot(t_vec, err_real, '-', 'Color', 'k', 'LineWidth', 1.2, ...
         'DisplayName', 'Real (With Meas. Error)');
    
    % Marcadores de Atualização
    sys_updates = sortrows(completed, 1);
    if ~isempty(sys_updates)
        plot(sys_updates(:,1), zeros(size(sys_updates,1),1), '^', ...
             'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerSize', 6, ...
             'DisplayName', 'Updates');
    end
    
    title('System Error Evolution: Real vs Ideal', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Position Error (m)', 'Interpreter', 'latex');
    legend('Location', 'best');
    
    % Limites (Zoom inicial)
    xlim([0, min(50, conf.Sim_Time)]); 
    ylim([0, max(err_real)*1.1]);
    
    % --- PLOT 2: Usuários Individuais (Apenas Real para não poluir) ---
    subplot(2, 1, 2); hold on; grid on; box on;
    title('Individual User Trajectories (Real Error)', 'Interpreter', 'latex', 'FontSize', 14);
    
    colors = lines(n_users_to_plot);
    for u = 1:n_users_to_plot
        mask = (completed(:, 3) == u);
        data_u = completed(mask, :);
        v_u = velocities(u);
        
        [t_u, err_u_real, ~] = generate_error_curve(data_u, v_u, conf.Sim_Time);
        plot(t_u, err_u_real, '-', 'Color', colors(u,:), 'LineWidth', 1.2, ...
             'DisplayName', sprintf('User %d', u));
    end
    
    xlabel('Time (s)', 'Interpreter', 'latex');
    ylabel('Error (m)', 'Interpreter', 'latex');
    xlim([0, min(50, conf.Sim_Time)]);
    legend('Location', 'best');
    
end

% =========================================================================
%  HELPERS
% =========================================================================

function [t_out, err_out_real, err_out_ideal] = generate_error_curve(completed_data, v, T_max)
    % Gera vetores de tempo e erro (Real e Ideal)
    
    data = sortrows(completed_data, 1);
    
    t_out = [];
    err_out_real = [];
    err_out_ideal = [];
    
    current_t = 0;
    last_gen = 0;
    last_err = 0; 
    
    dt_res = 0.05; % Resolução fina para ver o "V" shape
    
    finish_times = data(:, 1);
    gen_times    = data(:, 2);
    meas_errors  = data(:, 4); 
    
    for k = 1:length(finish_times)
        next_t = finish_times(k);
        if next_t > T_max, next_t = T_max; end
        
        if next_t > current_t
            t_seg = current_t : dt_res : next_t;
            if t_seg(end) ~= next_t, t_seg = [t_seg, next_t]; end
            
            % --- CÁLCULO CORE ---
            % 1. Ideal: O erro de medição é assumido como 0. 
            %    Erro = v * (t - t_gen). Sempre >= 0.
            %    (t_seg - last_gen) é a "Idade" (AoI).
            seg_ideal = abs(v * (t_seg - last_gen)); 
            
            % 2. Real: Inclui o last_err (pode ser negativo).
            %    Erro = |v * (t - t_gen) + erro_sensor|
            seg_real  = abs(v * (t_seg - last_gen) + last_err);
            
            t_out = [t_out, t_seg];
            err_out_ideal = [err_out_ideal, seg_ideal];
            err_out_real  = [err_out_real, seg_real];
        end
        
        current_t = next_t;
        last_gen = gen_times(k);
        last_err = meas_errors(k);
        
        if current_t >= T_max, break; end
    end
    
    % Resto até T_max
    if current_t < T_max
        t_seg = current_t : dt_res : T_max;
        
        seg_ideal = abs(v * (t_seg - last_gen));
        seg_real  = abs(v * (t_seg - last_gen) + last_err);
        
        t_out = [t_out, t_seg];
        err_out_ideal = [err_out_ideal, seg_ideal];
        err_out_real  = [err_out_real, seg_real];
    end
end

function S = get_erasure_service(n, k, delta, sym_dur)
    if delta >= 1, delta = 0.999; end 
    failures = nbinrnd(k, 1-delta, n, 1);
    S = (k + failures) * sym_dur;
end

function completed = run_lcfs_np_plot(data, T_max)
    % Recria a fila NP para obter os tempos de saída
    % data input: [Arr, Serv, Err, User] -> NÃO, o input no main é [Arr, User, Serv, Err]
    % Ajustando para ler colunas corretamente
    
    % Data structure esperada aqui baseada na chamada:
    % data_full = [all_arrivals(1,2), services(3), meas_errors(4)];
    % col 1: Time, col 2: User, col 3: Service, col 4: Error
    
    N = size(data, 1);
    completed = [];
    served_mask = false(N, 1);
    time_now = 0;
    
    while true
        if time_now >= T_max, break; end
        queue = find(data(:,1) <= time_now & ~served_mask);
        if isempty(queue)
            upcoming = find(data(:,1) > time_now & ~served_mask, 1);
            if isempty(upcoming), break; end
            time_now = data(upcoming, 1); continue;
        end
        target = queue(end);
        
        if length(queue) > 1, served_mask(queue(1:end-1)) = true; end
        
        fin_t = time_now + data(target, 3);
        
        if fin_t <= T_max
            % Output: [Finish, Gen, User, Error]
            completed = [completed; fin_t, data(target, 1), data(target, 2), data(target, 4)];
            time_now = fin_t;
        else
            time_now = fin_t; 
        end
        served_mask(target) = true;
    end
end