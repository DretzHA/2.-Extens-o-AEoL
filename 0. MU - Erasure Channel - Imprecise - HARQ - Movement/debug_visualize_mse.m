function debug_visualize_distance(updates, segments, T_max)
    % Função de Debug Visual: Decomposição do Erro com Áreas Coloridas e Pausa
    
    % --- Configuração da Figura ---
    h_fig = figure('Name', 'Decomposição do Erro: AoI + RMSE', 'Color', 'w', 'Position', [100 100 1000 800]);
    
    % Subplot 1: Mapa 2D
    ax1 = subplot(2,1,1); hold on; grid on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); title('Geometria do Erro: Movimento (AoI) + Sensor (RMSE)');
    % Trajetória de fundo
    plot(segments(:,3), segments(:,4), '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 2);
    
    % Subplot 2: Erro vs Tempo (Com Áreas)
    ax2 = subplot(2,1,2); hold on; grid on;
    xlabel('Tempo (s)'); ylabel('Erro (m)'); 
    title('Evolução do Erro: Área Vermelha = Erro Total | Área Azul = Componente de Atraso');
    xlim([0 T_max]);
    
    % --- Inicialização dos Objetos Gráficos (Handles) ---
    
    % -- Subplot 2 (Áreas e Linhas) --
    % Criamos os objetos de área vazios primeiro para ficarem no fundo
    h_area_total = area(ax2, 0, 0, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Área Erro Total');
    h_area_aoi   = area(ax2, 0, 0, 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Área Erro AoI');
    
    % Linhas de contorno
    h_line_total = plot(ax2, 0, 0, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Erro Total (AEoL)');
    h_line_aoi   = plot(ax2, 0, 0, 'b--', 'LineWidth', 1, 'DisplayName', 'Comp. Atraso');
    %legend(ax2, 'Location', 'northwest');
    
    % -- Subplot 1 (Mapa) --
    current_time = 0;
    p0 = get_position_at_t_debug(0, segments);
    estimated_pos = p0;
    pos_at_gen_truth = p0;
    
    h_real  = plot(ax1, p0(1), p0(2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Posição Real (t)');
    h_ghost = plot(ax1, p0(1), p0(2), 'o', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Posição na Geração (t_{gen})');
    h_est   = plot(ax1, estimated_pos(1), estimated_pos(2), 'rs', 'MarkerFaceColor', 'r', 'DisplayName', 'Estimativa Atual');
    
    h_vec_aoi   = plot(ax1, [p0(1) p0(1)], [p0(2) p0(2)], 'b--', 'DisplayName', 'Vetor AoI');
    h_vec_rmse  = plot(ax1, [p0(1) p0(1)], [p0(2) p0(2)], 'r-', 'LineWidth', 2, 'DisplayName', 'Vetor Sensor');
    h_vec_total = plot(ax1, [p0(1) p0(1)], [p0(2) p0(2)], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Vetor Total');
    
    % Buffers para acumular histórico do gráfico de tempo
    hist_t = [];
    hist_err_total = [];
    hist_err_aoi = [];
    
    % Ordenar eventos
    if ~isempty(updates)
        events = sortrows(updates, 1); 
    else
        events = [];
    end
    
    event_idx = 1;
    num_events = size(events, 1);
    
    % --- LOOP VISUAL ---
    while current_time < T_max
        if event_idx <= num_events
            next_event_time = min(events(event_idx, 1), T_max);
            is_update = true;
        else
            next_event_time = T_max;
            is_update = false;
        end
        
        if next_event_time > current_time
            t_start_step = current_time;
            
            while t_start_step < next_event_time
                % Lógica RWP para achar segmento
                idx = find(segments(:,1) <= t_start_step & segments(:,2) > t_start_step, 1, 'last');
                if isempty(idx), idx = find(segments(:,2) >= t_start_step, 1, 'first'); end
                if isempty(idx), break; end
                
                row = segments(idx, :);
                t_limit_step = min(next_event_time, row(2));
                
                % Vetor de tempo para suavidade
                t_vec = linspace(t_start_step, t_limit_step, 15); 
                
                chunk_total_err = zeros(length(t_vec), 1);
                chunk_aoi_err   = zeros(length(t_vec), 1);
                
                for k = 1:length(t_vec)
                    t_inst = t_vec(k);
                    
                    % 1. Onde o alvo está realmente (t atual)
                    dt = t_inst - row(1);
                    pos_real = [row(3) + row(5)*dt, row(4) + row(6)*dt];
                    
                    % 2. Erros
                    chunk_total_err(k) = norm(pos_real - estimated_pos);
                    chunk_aoi_err(k)   = norm(pos_real - pos_at_gen_truth);
                end
                
                % --- Atualiza Gráfico de Tempo (Acumula dados) ---
                hist_t = [hist_t, t_vec];
                hist_err_total = [hist_err_total; chunk_total_err];
                hist_err_aoi = [hist_err_aoi; chunk_aoi_err];
                
                % Atualiza objetos 'Area' e 'Line'
                set(h_area_total, 'XData', hist_t, 'YData', hist_err_total);
                set(h_area_aoi,   'XData', hist_t, 'YData', hist_err_aoi);
                set(h_line_total, 'XData', hist_t, 'YData', hist_err_total);
                set(h_line_aoi,   'XData', hist_t, 'YData', hist_err_aoi);
                
                % --- Atualiza Mapa 2D (Apenas o último ponto) ---
                set(h_real, 'XData', pos_real(1), 'YData', pos_real(2));
                set(h_ghost, 'XData', pos_at_gen_truth(1), 'YData', pos_at_gen_truth(2));
                
                % Vetores
                set(h_vec_aoi, 'XData', [pos_real(1) pos_at_gen_truth(1)], 'YData', [pos_real(2) pos_at_gen_truth(2)]);
                set(h_vec_rmse, 'XData', [pos_at_gen_truth(1) estimated_pos(1)], 'YData', [pos_at_gen_truth(2) estimated_pos(2)]);
                set(h_vec_total, 'XData', [pos_real(1) estimated_pos(1)], 'YData', [pos_real(2) estimated_pos(2)]);
                
                t_start_step = t_limit_step;
                drawnow; 
            end
        end
        
        current_time = next_event_time;
        
        % --- MOMENTO DA ATUALIZAÇÃO (PAUSA AQUI) ---
        if is_update && current_time < T_max
            gen_time = events(event_idx, 2);
            err_x = events(event_idx, 4);
            err_y = events(event_idx, 5);
            
            % 1. Visualiza a chegada do pacote (linha vertical)
            xline(ax2, current_time, 'k:', 'Alpha', 0.5);
            
            % Atualiza Título para avisar o usuário
            title(ax2, sprintf('UPDATE Recebido em t=%.2fs! (Gerado em t=%.2fs). PRESSIONE UMA TECLA...', current_time, gen_time), 'Color', 'r', 'FontWeight', 'bold');
            
            % --- PAUSA ---
            pause; 
            % -------------
            
            % Restaura título e processa
            title(ax2, 'Simulando... (Evolução do Erro)');
            
            % Atualiza a "Verdade Passada" e Estimativa
            pos_at_gen_truth = get_position_at_t_debug(gen_time, segments);
            new_estimated_pos = pos_at_gen_truth + [err_x, err_y];
            
            % Efeito visual instantâneo (Queda do erro) no gráfico de tempo
            % Adicionamos um ponto no mesmo tempo t, mas com o erro novo (para fazer o degrau vertical)
            new_err_total = norm(pos_at_gen_truth - new_estimated_pos); % Basicamente só o erro do sensor
            % O erro de AoI zera instantaneamente (pois t_now == t_gen para efeitos de tracking novo, se considerarmos update instantaneo, mas aqui t_gen é passado.
            % Na verdade, o AoI reseta para (current - gen), mas o "tracking" recomeça do ponto de geração.
            % Vamos visualizar o erro caindo para o erro do sensor.
            
            hist_t = [hist_t, current_time];
            hist_err_total = [hist_err_total; norm(get_position_at_t_debug(current_time, segments) - new_estimated_pos)];
            hist_err_aoi = [hist_err_aoi; norm(get_position_at_t_debug(current_time, segments) - pos_at_gen_truth)];
            
            estimated_pos = new_estimated_pos;
            set(h_est, 'XData', estimated_pos(1), 'YData', estimated_pos(2));
            
            event_idx = event_idx + 1;
        end
    end
    fprintf('Visualização concluída.\n');
end

function pos = get_position_at_t_debug(t, segments)
    idx = find(segments(:,1) <= t & segments(:,2) >= t, 1, 'last');
    if isempty(idx)
        idx = size(segments, 1); 
        t = min(t, segments(idx, 2));
    end
    row = segments(idx, :);
    dt = t - row(1);
    pos = [row(3) + row(5)*dt, row(4) + row(6)*dt];
end