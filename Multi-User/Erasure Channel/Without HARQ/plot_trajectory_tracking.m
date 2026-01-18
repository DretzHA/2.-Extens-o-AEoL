function plot_trajectory_tracking_comparison(conf)
    % Gera comparação visual de Rastreamento de Trajetória: NP vs Preemptivo
    % Baseado na mesma semente de chegadas/erros para comparação justa.
    fprintf('Gerando comparação de Rastreamento (NP vs P)...\n');
    
    % --- 1. Configuração e Geração de Dados (Comum aos dois) ---
    T_plot = 50; % Janela de tempo para visualização
    
    % Usar a maior carga para garantir bastante atividade
    % Se der erro de índice, mude para 'end' ou um número menor
    idx_lambda = 5; 
    lambda_tot = conf.lambda_total_vec(idx_lambda); 
    lambda_user = lambda_tot / conf.num_users;
    
    % Estimar numero de pacotes
    est_n_pkts = ceil(T_plot * lambda_user * 2) + 20;
    
    % Gerar Chegadas
    dt = exprnd(1/lambda_user, est_n_pkts, 1);
    t_arr = cumsum(dt);
    t_arr = t_arr(t_arr <= T_plot);
    
    if isempty(t_arr), warning('Nenhum pacote gerado.'); return; end
    n_pkts = length(t_arr);
    
    % Gerar Serviços e Erros (Idênticos para NP e P)
    fixed_service = conf.k * conf.Symbol_Duration;
    services = ones(n_pkts, 1) * fixed_service;
    pkt_success_prob = (1 - conf.delta)^conf.k;
    is_successful = rand(n_pkts, 1) < pkt_success_prob;
    
    meas_errors = normrnd(0, conf.sensor_sigma, n_pkts, 1);
    
    v = conf.velocities(1); 
    
    % Estrutura de dados: [Arr, User(dummy), Serv, Err]
    data_in = [t_arr, ones(n_pkts, 1), services, meas_errors, is_successful];
    
    % --- 2. Rodar as Filas ---
    completed_np = run_lcfs_np_tracking(data_in, T_plot);
    completed_p  = run_lcfs_p_tracking(data_in, T_plot);
    
    % --- 3. Plotagem ---
    fig = figure('Name', 'Trajectory Tracking: NP vs P', 'Color', 'w', 'Position', [100 50 1000 800]);
    
    % >> SUBPLOT 1: NP <<
    subplot(2, 1, 1); hold on; grid on; box on;
    plot_single_tracking_subplot(completed_np, v, T_plot, 'LCFS - Non Preemptive (NP)');
    
    % >> SUBPLOT 2: P <<
    subplot(2, 1, 2); hold on; grid on; box on;
    plot_single_tracking_subplot(completed_p, v, T_plot, 'LCFS - Preemptive (P)');
    
    % Hack para legenda não duplicada (opcional)
    h = axes(fig, 'visible','off'); 
end

function plot_single_tracking_subplot(completed_data, v, T_max, title_str)
    % Função auxiliar atualizada com marcadores de GERAÇÃO
    
    if isempty(completed_data)
        text(T_max/2, v*T_max/2, 'No Updates Processed', 'HorizontalAlignment', 'center');
        title(title_str, 'Interpreter', 'latex', 'FontSize', 14);
        return;
    end
    completed_data = sortrows(completed_data, 1);
    
    % Preparar vetores para 'stairs'
    update_times = [0; completed_data(:, 1)];
    gen_times    = [0; completed_data(:, 2)];
    errors       = [0; completed_data(:, 4)];
    
    % Cálculo das estimativas
    val_ideal = v * gen_times;           
    val_real  = v * gen_times + errors;  
    
    % Fechar o gráfico até T_max
    update_times = [update_times; T_max];
    val_ideal    = [val_ideal; val_ideal(end)];
    val_real     = [val_real; val_real(end)];
    
    % --- Plotting Curves ---
    % 1. Ground Truth (Linha Sólida Vinho)
    fplot(@(t) v*t, [0 T_max], 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.5, ...
        'DisplayName', 'Ground Truth');
    
    % 2. Estimativa Ideal (Pontilhada Verde Clara)
    stairs(update_times, val_ideal, ':', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2.0, ...
        'DisplayName', 'Precise Est.');
    
    % 3. Estimativa Real (Sólida Verde Escura)
    stairs(update_times, val_real, '-', 'Color', [0 0.5 0], 'LineWidth', 2.0, ...
        'DisplayName', 'Imprecise Est.');
    
    % 4. Area de Erro (Visual)
    t_fine = linspace(0, T_max, 500);
    pos_true = v * t_fine;
    pos_est  = interp1(update_times, val_real, t_fine, 'previous');
    fill([t_fine, fliplr(t_fine)], [pos_true, fliplr(pos_est)], ...
         [0 0.5 0], 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
     
    % --- MARCADORES DE GERAÇÃO E ATUALIZAÇÃO ---
    
    % Selecionar apenas alguns pontos para não poluir o gráfico
    % Mostra os primeiros 15 ou menos
    idx_show = 1:min(3, length(completed_data));
    
    if ~isempty(completed_data)
        data_show = completed_data(idx_show, :);
        
        % A. Marcadores de GERAÇÃO (Triângulo Vermelho para cima)
        % Eixo X: Tempo de Geração (Col 2)
        % Eixo Y: Posição Real na Geração (v * Tempo de Geração)
        t_gen = data_show(:, 2);
        p_gen = v * t_gen;
        
        plot(t_gen, p_gen, '^', ...
            'MarkerFaceColor', [0.8500 0.3250 0.0980], ... % Vermelho/Laranja
            'MarkerEdgeColor', 'none', 'MarkerSize', 5, ...
            'DisplayName', 'Gen. Point');
            
        % B. Marcadores de ATUALIZAÇÃO (Triângulo Verde para baixo)
        % Eixo X: Tempo de Chegada (Col 1)
        % Eixo Y: Valor da Estimativa (v * Gen + Erro)
        t_upd = data_show(:, 1);
        p_upd = v * t_gen + data_show(:, 4); % Nota: usa t_gen para calcular P, mas plota em t_upd
        
        plot(t_upd, p_upd, 'v', ...
            'MarkerFaceColor', [0 0.5 0], ... % Verde Escuro
            'MarkerEdgeColor', 'none', 'MarkerSize', 5, ...
            'DisplayName', 'Update Point');
            
        % C. Linha de Conexão (Delay Visual)
        % Conecta o ponto de Geração ao ponto de Atualização correspondente
        for k = 1:length(t_gen)
            plot([t_gen(k) t_upd(k)], [p_gen(k) p_upd(k)], ...
                ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end

    % Formatação
    title(title_str, 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Position (m)', 'Interpreter', 'latex');
    xlabel('Time (s)', 'Interpreter', 'latex');
    xlim([0 T_max]);
    ylim([0, v*T_max * 1.1]); 
    grid minor;
    legend('Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
end

% =========================================================================
%  FILAS (Tracking Logic) - MANTIDAS IGUAIS
% =========================================================================
function completed = run_lcfs_np_tracking(data, T_max)
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
            % Check Success
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


function completed = run_lcfs_p_tracking(data, T_max)
    N = size(data, 1); completed = [];
    for k = 1:N-1
        arr_c = data(k,1); arr_n = data(k+1,1); serv = data(k,3);
        finish_time = arr_c + serv;
        
        if finish_time <= arr_n && finish_time <= T_max
             % Check Success
             if data(k, 5) == 1
                 completed = [completed; finish_time, arr_c, data(k,2), data(k,4)];
             end
        end
    end
    % Last pkt
    finish_time = data(N,1) + data(N,3);
    if finish_time <= T_max
        if data(N, 5) == 1
            completed = [completed; finish_time, data(N,1), data(N,2), data(N,4)];
        end
    end
end
