function debug_visualize_mse(updates, segments, T_max)
    % Função de Debug Visual: Decomposição do Erro + AoI + Marcadores de Geração
    
    h_fig = figure('Name', 'Validação Completa: Erro + AoI', 'Color', 'w', 'Position', [100 50 1000 900]);
    
    % === SUBPLOT 1: Mapa 2D ===
    ax1 = subplot(3,1,1); hold on; grid on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); title('1. Mapa');
    plot(segments(:,3), segments(:,4), '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 2);
    xlim([0 60]); ylim([0 60]); % Descomente se quiser fixar o zoom
    
    % === SUBPLOT 2: Erro (AEoL) ===
    ax2 = subplot(3,1,2); hold on; grid on;
    xlabel('Tempo (s)'); ylabel('Erro (m)'); 
    title('2. AEoL');
    xlim([0 T_max]);
    
    % === SUBPLOT 3: Age of Information (AoI) ===
    ax3 = subplot(3,1,3); hold on; grid on;
    xlabel('Tempo (s)'); ylabel('AoI (s)');
    title('3. AoI ');
    xlim([0 T_max]);
    
    % --- Inicialização dos Objetos Gráficos ---
    
    % Subplot 2
    h_area_total = area(ax2, 0, 0, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h_area_aoi   = area(ax2, 0, 0, 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h_line_total = plot(ax2, 0, 0, 'r-', 'LineWidth', 1.5);

    legend(ax2, [h_area_total, h_area_aoi], ...
           'Erro Total', ...       % Legenda do Vermelho
           'Erro por AoI', ... % Legenda do Azul
           'Location', 'northwest', 'AutoUpdate', 'off');
    
    % Subplot 3
    h_area_val_aoi = area(ax3, 0, 0, 'FaceColor', [0.2 0.8 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'k');
    % Linha base no zero para referência dos marcadores
    yline(ax3, 0, 'k-', 'Alpha', 0.2); 
    
    % Subplot 1
    p0 = get_position_at_t_debug(0, segments);
    estimated_pos = p0;
    pos_at_gen_truth = p0;
    
    h_real  = plot(ax1, p0(1), p0(2), 'bo', 'MarkerFaceColor', 'b');
    h_est   = plot(ax1, estimated_pos(1), estimated_pos(2), 'rs', 'MarkerFaceColor', 'r');
    h_vec_total = plot(ax1, [p0(1) p0(1)], [p0(2) p0(2)], 'k:', 'LineWidth', 1.5);
    
    % Variáveis de Estado
    current_time = 0;
    last_gen_time = 0; 
    
    hist_t = [];
    hist_err_total = [];
    hist_err_aoi = [];
    hist_val_aoi = [];
    
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
                % Achar segmento
                idx = find(segments(:,1) <= t_start_step & segments(:,2) > t_start_step, 1, 'last');
                if isempty(idx), idx = find(segments(:,2) >= t_start_step, 1, 'first'); end
                if isempty(idx), break; end
                
                row = segments(idx, :);
                t_limit_step = min(next_event_time, row(2));
                
                t_vec = linspace(t_start_step, t_limit_step, 15); 
                
                chunk_total_err = zeros(length(t_vec), 1);
                chunk_err_aoi_comp = zeros(length(t_vec), 1);
                chunk_val_aoi   = zeros(length(t_vec), 1);
                
                for k = 1:length(t_vec)
                    t_inst = t_vec(k);
                    dt = t_inst - row(1);
                    pos_real = [row(3) + row(5)*dt, row(4) + row(6)*dt];
                    
                    chunk_total_err(k) = norm(pos_real - estimated_pos);
                    chunk_err_aoi_comp(k) = norm(pos_real - pos_at_gen_truth);
                    chunk_val_aoi(k) = t_inst - last_gen_time;
                end
                
                % Atualiza Arrays
                hist_t = [hist_t, t_vec];
                hist_err_total = [hist_err_total; chunk_total_err];
                hist_err_aoi = [hist_err_aoi; chunk_err_aoi_comp];
                hist_val_aoi = [hist_val_aoi; chunk_val_aoi];
                
                % Update Plots
                set(h_area_total, 'XData', hist_t, 'YData', hist_err_total);
                set(h_area_aoi,   'XData', hist_t, 'YData', hist_err_aoi);
                set(h_line_total, 'XData', hist_t, 'YData', hist_err_total);
                set(h_area_val_aoi, 'XData', hist_t, 'YData', hist_val_aoi);
                
                set(h_real, 'XData', pos_real(1), 'YData', pos_real(2));
                set(h_est, 'XData', estimated_pos(1), 'YData', estimated_pos(2));
                set(h_vec_total, 'XData', [pos_real(1) estimated_pos(1)], 'YData', [pos_real(2) estimated_pos(2)]);
                
                t_start_step = t_limit_step;
                drawnow; 
            end
        end
        
        current_time = next_event_time;
        
        % --- UPDATE CHEGOU ---
        if is_update && current_time < T_max
            gen_time = events(event_idx, 2);
            err_x = events(event_idx, 4);
            err_y = events(event_idx, 5);
            
            % Linhas verticais de evento
            xline(ax2, current_time, 'k:', 'Alpha', 0.5);
            xline(ax3, current_time, 'k:', 'Alpha', 0.5);
            
            % === NOVO: PLOTAR MARCADOR DE T_GEN NO EIXO DO AOI ===
            % Triângulo verde no eixo X indicando quando nasceu
            plot(ax3, gen_time, 0, '^', 'MarkerEdgeColor', 'k', ...
                 'MarkerFaceColor', [0 0.7 0], 'MarkerSize', 8, 'DisplayName', 't_{gen}');
             
            % Opcional: Texto indicando o tempo
            text(ax3, gen_time, 0, sprintf(' %.1fs', gen_time), ...
                 'VerticalAlignment', 'top', 'FontSize', 7, 'Color', [0 0.5 0]);
            
            % Feedback Visual
            % title(ax2, sprintf('UPDATE! Chegou em %.2fs (Gerado em %.2fs)', current_time, gen_time), 'Color', 'r');
            pause;
            % title(ax2, 'Evolução do Erro (AEoL)');
            
            % Cálculos de Posição
            pos_at_gen_truth = get_position_at_t_debug(gen_time, segments);
            new_estimated_pos = pos_at_gen_truth + [err_x, err_y];
            
            last_gen_time = gen_time;
            
            % Degrau visual nos gráficos
            hist_t = [hist_t, current_time];
            
            % O Erro muda
            hist_err_total = [hist_err_total; norm(get_position_at_t_debug(current_time, segments) - new_estimated_pos)];
            hist_err_aoi = [hist_err_aoi; norm(get_position_at_t_debug(current_time, segments) - pos_at_gen_truth)];
            
            % O AoI cai para o "System Delay" (Current - Gen)
            new_val_aoi = current_time - last_gen_time;
            hist_val_aoi = [hist_val_aoi; new_val_aoi];
            
            % === NOVO: LINHA DE CONEXÃO VISUAL (DELAY) ===
            % Desenha uma linha pontilhada verde do t_gen até o novo valor do AoI
            plot(ax3, [gen_time, current_time], [0, new_val_aoi], 'g:', 'LineWidth', 1);
            
            estimated_pos = new_estimated_pos;
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