function plot_sawtooth_woHARQ(conf, n_users_to_plot)
    % Gera uma realização única e plota Trajetória Real vs Ideal
    % NO HARQ Version
    
    fprintf('Gerando gráfico Sawtooth (Erasure Only)...\\n');
    
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
    
    % --- MUDANÇA: Serviço Fixo e Bernoulli Sucesso ---
    fixed_service = conf.k * conf.Symbol_Duration;
    services = ones(n_pkts, 1) * fixed_service;
    
    pkt_success_prob = (1 - conf.delta)^conf.k;
    is_successful = rand(n_pkts, 1) < pkt_success_prob;
    
    meas_errors = normrnd(0, conf.sensor_sigma, n_pkts, 1);
    
    % Data structure: [Arr, User, Serv, Err, Success]
    data_plot = [all_arrivals, services, meas_errors, is_successful];
    
    % Rodar NP só para pegar os pontos de update
    completed = run_lcfs_np_plot(data_plot, conf.Sim_Time);
    
    % --- Plotting ---
    figure('Name', 'Sawtooth Trajectory (Erasure)', 'Color', 'w', 'Position', [100 100 1000 400]);
    
    % Plotar apenas para o usuário 1 para clareza
    u_plot = 1;
    updates_u = completed(completed(:,3) == u_plot, :);
    
    % Trajetória Real (Linha contínua crescendo e caindo)
    % Vamos reconstruir a curva de erro ponto a ponto
    t_grid = 0:0.1:conf.Sim_Time;
    err_traj = zeros(size(t_grid));
    
    last_gen = 0; 
    last_err = 0;
    idx_up = 1;
    v = conf.velocities(u_plot);
    
    for i = 1:length(t_grid)
        t = t_grid(i);
        % Verifica se houve update
        if idx_up <= size(updates_u, 1)
            if t >= updates_u(idx_up, 1)
                last_gen = updates_u(idx_up, 2);
                last_err = updates_u(idx_up, 4);
                idx_up = idx_up + 1;
            end
        end
        % Erro atual: | v*(t - t_gen) + noise |
        err_traj(i) = abs(v * (t - last_gen) + last_err);
    end
    
    plot(t_grid, err_traj, 'b-', 'LineWidth', 1.5); hold on;
    
    % Marcar Updates
    if ~isempty(updates_u)
        plot(updates_u(:,1), zeros(size(updates_u,1),1), 'g^', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
    end
    
    xlabel('Time (s)'); ylabel('Position Error (m)');
    title(sprintf('Error Trajectory User 1 (Erasure Channel, P_{succ}=%.2f)', pkt_success_prob));
    grid on;
    xlim([0 100]); % Zoom no início
end

function completed = run_lcfs_np_plot(data, T_max)
    % Versão simplificada local para o plot
    N = size(data, 1); completed = []; served_mask = false(N, 1); time_now = 0;
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
            % Só salva se sucesso
            if data(target, 5) == 1
                completed = [completed; fin_t, data(target, 1), data(target, 2), data(target, 4)];
            end
            time_now = fin_t;
        else
            time_now = fin_t; 
        end
        served_mask(target) = true;
    end
end
